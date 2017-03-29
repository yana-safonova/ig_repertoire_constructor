#!/usr/bin/env python2

import os
import sys
import init
import logging
import shutil
import ntpath

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) #+ '/'
py_src = os.path.join(home_directory, "py/pipeline/")
vjf_config_dir = os.path.join(home_directory, "configs/vj_finder")
vj_finder_bin = os.path.join(home_directory, "build/release/bin/vj_finder")
run_vj_finder = os.path.join(home_directory, "build/release/bin/./vj_finder")

sys.path.append(py_src)
import process_cfg
import support

test_reads = os.path.join(home_directory, "test_dataset/merged_reads.fastq")
test_dir = os.path.join(home_directory, "vjf_test")

tool_name = "VJ Finder"

################################################################################
def CheckBinariesExistance(log):
    if not os.path.exists(vj_finder_bin):
        log.info("ERROR: " + tool_name + " binary was not found. Please type make before running")
        sys.exit(1)

################################################################################
def CheckInputReads(input_reads, log):
    if not os.path.exists(input_reads):
        log.info("ERROR: input reads " + input_reads + " were not found")
        sys.exit(1)

def CheckGermlineLocus(locus, log):
    loci_set = ['IGH', 'IGK', 'IGL', 'IG', 'TRA', 'TRB', 'TRG', 'TRD', "TR", "all"]
    if locus not in loci_set:
        log.info("ERROR: locus option (-l) " + locus + " was not recognized. "
                                                       "Please use one of the available option values: " + str(loci_set))
        sys.exit(1)

def CheckOrganism(organism, log):
    org_set = ['human', 'mouse', 'rat', 'pig', 'rabbit', 'rhesus_monkey']
    if organism not in org_set:
        log.info("ERROR: organism option " + organism + " was not recognized. "
                                                        "Please use one of the available option values: " + str(org_set))
        sys.exit(1)

def CheckParamsCorrectness(params, log):
    CheckInputReads(params.input_reads, log)
    CheckGermlineLocus(params.loci, log)
    CheckOrganism(params.organism, log)

################################################################################
def SetOutputParams(params, log):
    if params.input_reads == test_reads:
        params.output_dir = test_dir
    if params.input_reads != test_reads and params.output_dir == "":
        log.info("ERROR: Output dir (-o) was not specified")
        sys.exit(1)
    params.output_dir = os.path.abspath(params.output_dir)
    params.config_dir = os.path.join(params.output_dir, "configs")
    params.config_file = os.path.join(vjf_config_dir, "config.info")

def PrepareOutputDir(params):
    if os.path.exists(params.output_dir):
        shutil.rmtree(params.output_dir)
    os.mkdir(params.output_dir)

def PrintParams(params, log):
    log.info(tool_name + " parameters:")
    log.info("  Input reads:\t\t" + params.input_reads)
    log.info("  Output directory:\t" + params.output_dir)
    log.info("  Numner of threads:\t" + str(params.num_threads) + "\n")

    log.info("  Loci:\t\t\t" + params.loci)
    log.info("  Organism:\t\t" + params.organism)
    log.info("  Using pseudogenes:\t" + str(not params.no_pseudogenes) + "\n")

    log.info("  Size of k-mer for V alignment:\t\t" + str(params.v_kmer))
    log.info("  Size of k-mer for J alignment:\t\t" + str(params.j_kmer) + "\n")

################################################################################
def main(argv):
    from argparse import ArgumentParser
    parser = ArgumentParser(description="== " + tool_name + ": a tool for VJ alignment of full-length Rep-seq reads ==",
                            epilog="In case you have troubles running " + tool_name + ", you can write to igtools_support@googlegroups.com."
                            "Please provide us with vj_finder.log file from the output directory.",
                            add_help=False)
    req_args = parser.add_argument_group("Required params")
    input_args = req_args.add_mutually_exclusive_group(required=True)
    input_args.add_argument("-i", "--input",
                            type=str,
                            default="",
                            dest="input_reads",
                            help="Input reads in FASTQ/FASTA format")
    input_args.add_argument("--test",
                            action="store_const",
                            const=test_reads,
                            dest="input_reads",
                            help="Running test dataset")

    out_args = parser.add_argument_group("Output")
    out_args.add_argument("-o", "--output",
                          type=str,
                          dest="output_dir",
                          default="", #os.path.join(home_directory, "cdr_test"),
                          help="Output directory")

    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument("-t", "--threads",
                               type=int,
                               default=16,
                               dest="num_threads",
                               help="Threads number [default: %(default)d]")
    optional_args.add_argument("-h", "--help",
                               action="help",
                               help="Help message and exit")

    germline_args = parser.add_argument_group("Germline arguments")
    germline_args.add_argument("-l", "--loci",
                               type=str,
                               default="all",
                               dest="loci",
                               help="Loci: IGH, IGK, IGL, IG (all BCRs), TRA, TRB, TRG, TRD, TR (all TCRs) or all. "
                                    "[default: %(default)s]")
    germline_args.add_argument("-s", "--organism",
                               type=str,
                               default="human",
                               dest="organism",
                               help="Organism: human, mouse, pig, rabbit, rat, rhesus_monkey are available. "
                               "[default: %(default)s]")
    germline_args.add_argument("--no-pseudogenes",
                                action='store_const',
                                const=True,
                                dest="no_pseudogenes",
                                help = "Exclusion of pseudogenes for alignment of input reads")

    alignment_args = parser.add_argument_group("Alignment arguments")
    alignment_args.add_argument("--v-kmer",
                                type=int,
                                default=7,
                                dest="v_kmer",
                                help="Size of k-mer that is used for V gene alignment")
    alignment_args.add_argument("--j-kmer",
                                type=int,
                                default=5,
                                dest="j_kmer",
                                help="Size of k-mer that is used for J gene alignment")
    alignment_args.add_argument("--v-min-cov",
                                type=int,
                                default=50,
                                dest="v_min_coverage",
                                help="Length of minimal coverage of V gene by k-mers")
    alignment_args.add_argument("--j-min-cov",
                                type=int,
                                default=13,
                                dest="j_min_coverage",
                                help="Length of minimal coverage of J gene by k-mers")
    #alignment_args.add_argument("--v-num-cand",
    #                            type=int,
    #                            default=3,
    #                            dest="v_min_candidates",
    #                            help="Number of candidates reported ")

    # prepare log
    log = logging.getLogger('vj_finder')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    params = parser.parse_args()

    CheckBinariesExistance(log)
    CheckParamsCorrectness(params, log)
    SetOutputParams(params, log)
    PrepareOutputDir(params)

    # log file
    params.log_filename = os.path.join(params.output_dir, "diversity_analyzer.log")
    if os.path.exists(params.log_filename):
         os.remove(params.log_filename)
    log_handler = logging.FileHandler(params.log_filename, mode='a')
    log.addHandler(log_handler)

    # print command line
    command_line = "Command_line: "
    command_line += " ".join(argv)
    log.info(command_line + "\n")
    PrintParams(params, log)
    log.info("Log will be written to " + params.log_filename + "\n")

    # PrepareConfigs(params, log)
    # try:
    #     cdr_command_line = run_cdr_labeler + " " + params.cdr_labeler_config_file
    #     support.sys_call(cdr_command_line, log)
    #     if not params.skip_plots:
    #         log.info("\n==== Visualization of diversity statistics ====")
    #         visualize_vj_stats.main(["", os.path.join(params.output_dir, "cdr_details.txt"),
    #                              os.path.join(params.output_dir, "shm_details.txt"),
    #                              os.path.join(params.output_dir, "plots"), log])
    #         log.info("\n==== Annotation report creation ====")
    #         html_report_writer.main(os.path.join(params.output_dir, "cdr_details.txt"),
    #                             os.path.join(params.output_dir, "shm_details.txt"),
    #                             os.path.join(params.output_dir, "plots"),
    #                             os.path.join(params.output_dir, "annotation_report.html"), log)
    #     Cleanup(params, log)
    #     log.info("\nThank you for using " + tool_name + "!\n")
    # except (KeyboardInterrupt):
    #      log.info("\n" + tool_name + " was interrupted!")
    # except Exception:
    #     exc_type, exc_value, _ = sys.exc_info()
    #     if exc_type == SystemExit:
    #         sys.exit(exc_value)
    #     else:
    #         log.exception(exc_value)
    #         log.info("\nERROR: Exception caught.")
    #          #supportInfo(log)
    # except BaseException:
    #     exc_type, exc_value, _ = sys.exc_info()
    #     if exc_type == SystemExit:
    #         sys.exit(exc_value)
    #     else:
    #         log.exception(exc_value)
    #         log.info("\nERROR: Exception caught.")
    #         #supportInfo(log)
    #
    # log.info("Log was written to " + params.log_filename)


if __name__ == '__main__':
    main(sys.argv)