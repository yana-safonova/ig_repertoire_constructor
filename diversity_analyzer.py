#!/usr/bin/env python2

import os
import sys
import init
import logging
import shutil
import ntpath

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'
py_src = os.path.join(home_directory, "py/pipeline/")
cdr_labeler_config_dir = os.path.join(home_directory, "configs/cdr_labeler/")
vj_finder_config_dir = os.path.join(home_directory, "configs/vj_finder")
cdr_labeler_bin = os.path.join(home_directory, "build/release/bin/cdr_labeler")
hgc_bin = os.path.join(home_directory, "build/release/bin/ig_swgraph_construct")
run_cdr_labeler = os.path.join(home_directory, "build/release/bin/./cdr_labeler")
visualizer_dir = os.path.join(home_directory, "py/diversity_stats_visualizer")
data_annotation_dir = os.path.join(home_directory, "data/annotation")

sys.path.append(py_src)
sys.path.append(visualizer_dir)

import process_cfg
import support

import visualize_vj_stats
import html_report_writer

test_reads = os.path.join(home_directory, "test_dataset/merged_reads.fastq")
test_dir = "divan_test"

tool_name = "Diversity Analyzer"

def CheckBinariesExistance(params, log):
    if not os.path.exists(cdr_labeler_bin) or not os.path.exists(hgc_bin):
        log.info("ERROR: Binary files were not found. Please compile " + tool_name + " before running.")
        sys.exit(1)

def DomainParamCorrect(domain_str):
    return domain_str == "imgt" or domain_str == "kabat"

def LociParamIsIg(loci_str):
    return loci_str == "IG" or loci_str == "IGH" or loci_str == "IGK" or loci_str == "IGL"

organism_dict = {'human' : 'human', 'mouse' : 'mouse', 'rat' : 'rat',
                 'rabbit' : 'rabbit', 'rhesus-monkey' : 'rhesus_monkey'}

def OrganismParamCorrect(org_str):
    return org_str in organism_dict

def CheckParamsCorrectness(params, log):
    if not os.path.exists(params.input_reads):
        log.info("Input reads " + params.input_reads + " were not found")
        sys.exit(1)
    if not DomainParamCorrect(params.domain_system):
        log.info("Domain system " + params.domain_system + " is not recognized")
        sys.exit(1)
#    if not LociParamCorrect(params.loci):
#        log.info("Loci " + params.loci + " is not recognized")
#        sys.exit(1)
    if not OrganismParamCorrect(params.organism):
        log.info("Organism " + params.organism + " is not recognized")
        sys.exit(1)
    else:
        params.organism = organism_dict[params.organism]

def SetOutputParams(params, log):
    if params.input_reads == test_reads:
        params.output_dir = test_dir
    if params.input_reads != test_reads and params.output_dir == "":
        log.info("ERROR: Output dir (-o) was not specified")
        sys.exit(1)
    if not LociParamIsIg(params.loci):
        params.skip_plots = True
    params.output_dir = os.path.abspath(params.output_dir)
    params.config_dir = os.path.join(params.output_dir, "configs")
    params.cdr_config_file = os.path.join(cdr_labeler_config_dir, "config.info")
    params.vj_finder_config_file = os.path.join(vj_finder_config_dir, "config.info")

def PrepareOutputDir(params):
    if os.path.exists(params.output_dir):
        shutil.rmtree(params.output_dir)
    os.mkdir(params.output_dir)

def PrintParams(params, log):
    log.info(tool_name + " parameters:")
    log.info("  Input reads:\t\t" + params.input_reads)
    log.info("  Output directory:\t" + params.output_dir + "\n")
    log.info("  Domain system:\t" + params.domain_system)
    log.info("  Loci:\t\t\t" + params.loci)
    log.info("  Organism:\t\t" + params.organism + "\n")

########################################################################################################################

def CopyConfigs(params, log):
    if os.path.exists(params.config_dir):
        shutil.rmtree(params.config_dir)
    os.mkdir(params.config_dir)
    params.cdr_labeler_config_dir = os.path.abspath(os.path.join(params.config_dir, "cdr_labeler"))
    shutil.copytree(cdr_labeler_config_dir, params.cdr_labeler_config_dir)
    params.vj_finder_config_dir = os.path.abspath(os.path.join(params.config_dir, "vj_finder"))
    shutil.copytree(vj_finder_config_dir, params.vj_finder_config_dir)
    params.vj_finder_config_file = os.path.join(params.vj_finder_config_dir, "config.info")
    if not os.path.exists(params.vj_finder_config_file):
        log.info("ERROR: Config file " + params.vj_finder_config_file + " was not found")
        sys.exit(1)
    params.cdr_labeler_config_file = os.path.join(params.cdr_labeler_config_dir, "config.info")
    if not os.path.exists(params.cdr_labeler_config_file):
        log.info("ERROR: Config file " + params.cdr_labeler_config_file + " was not found")
        sys.exit(1)

def UpdateGermlineConfigFile(params, log):
    lines = open(params.germline_config_file, "r").readlines()
    fhandler = open(params.germline_config_file, "w")
    fhandler.write(lines[0].strip() + "\n")
    for i in range(1, len(lines)):
        splits = lines[i].strip().split()
        fhandler.write(splits[0] + "\t" + splits[1] + "\t" + os.path.join(home_directory, splits[2]) + "\n")
    fhandler.close()

def ModifyParamsWrtOrganism(params, cdr_param_dict, vj_param_dict):
    cdr_param_dict['imgt_v_annotation'] = os.path.join(data_annotation_dir, params.organism + "_v_imgt.txt")
    cdr_param_dict['kabat_v_annotation'] = os.path.join(data_annotation_dir, params.organism + "_v_kabat.txt")
    cdr_param_dict['imgt_j_annotation'] = os.path.join(data_annotation_dir, params.organism + "_j_imgt.txt")
    cdr_param_dict['kabat_j_annotation'] = os.path.join(data_annotation_dir, params.organism + "_j_kabat.txt")
    vj_param_dict['organism'] = params.organism

def ModifyConfigFiles(params, log):
    cdr_param_dict = dict()
    cdr_param_dict['input_reads'] = params.input_reads
    cdr_param_dict['output_dir'] = params.output_dir
    cdr_param_dict['vj_finder_config'] = params.vj_finder_config_file
    cdr_param_dict['num_threads'] = params.num_threads
    cdr_param_dict['domain_system'] = params.domain_system
    cdr_param_dict['run_hg_constructor'] = os.path.join(home_directory, './build/release/bin/ig_swgraph_construct')

    vj_param_dict = dict()
    vj_param_dict['loci'] = params.loci
    vj_param_dict['germline_dir'] = os.path.join(home_directory, "data/germline")
    params.germline_config_file = os.path.join(params.vj_finder_config_dir, "germline_files_config.txt")
    vj_param_dict['germline_filenames_config'] = params.germline_config_file

    ModifyParamsWrtOrganism(params, cdr_param_dict, vj_param_dict)
    process_cfg.substitute_params(params.cdr_labeler_config_file, cdr_param_dict, log)
    process_cfg.substitute_params(params.vj_finder_config_file, vj_param_dict, log)

def PrepareConfigs(params, log):
    CopyConfigs(params, log)
    ModifyConfigFiles(params, log)
    #UpdateGermlineConfigFile(params, log)

def Cleanup(params, log):
    params.trash_log = os.path.join(params.output_dir, "trash.out")
    if os.path.exists(params.trash_log):
        os.remove(params.trash_log)
    params.cdr3_graph = os.path.join(params.output_dir, "cdr3_graph.graph")
    if os.path.exists(params.cdr3_graph):
        os.remove(params.cdr3_graph)

########################################################################################################################

def main(argv):
    from argparse import ArgumentParser
    parser = ArgumentParser(description="== " + tool_name + ": a tool for diversity analysis of full-length immunosequencing reads ==",
                            epilog="In case you have troubles running " + tool_name + ", you can write to igtools_support@googlegroups.com."
                            "Please provide us with diversity_analyzer.log file from the output directory.",
                            add_help=False)
    req_args = parser.add_argument_group("Required params")
    input_args = req_args.add_mutually_exclusive_group(required=True)
    input_args.add_argument("-i", "--input",
                            type=str,
                            default="",
                            dest="input_reads",
                            help="Input reads in FASTQ/FATSA format")
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
                               help="Threads number. [default: %(default)d]")
    optional_args.add_argument("-d", "--domain",
                               type=str,
                               default="imgt",
                               dest="domain_system",
                               help="Domain system for CDR search: imgt OR kabat. [default: %(default)s]")

    vj_finder_args= parser.add_argument_group("VJ alignment params")
    optional_args.add_argument("-l", "--loci",
                               type=str,
                               default="all",
                               dest="loci",
                               help="Loci: IGH, IGK, IGL, IG (all BCRs), TRA, TRB, TRG, TRD, TR (all TCRs) or all. "
                                    "[default: %(default)s]")

    optional_args.add_argument("--org",
                               type=str,
                               default="human",
                               dest="organism",
                               help="Organism: human, mouse, rat, rabbit, rhesus-monkey. [default: %(default)s]")

    optional_args.add_argument("--skip-plots",
                               action="store_const",
                               const=True,
                               dest = "skip_plots",
                               help = "Skip drawing plots")

    optional_args.add_argument("--preset",
                               type=str,
                               default="divan",
                               dest="preset",
                               help="Predefined set of repertoire features to report. Supported values: divan."
                                    # "min outputs same fields as MiXCR export with min preset, "
                                    "divan outputs default IgDiversityAnalyzer feature set, "
                                    # "custom reads feature set from cdr_labeler config file."
                                    " [default: %(default)s]")

    optional_args.add_argument("-h", "--help",
                               action="help",
                               help="Help message and exit")

    # prepare log
    log = logging.getLogger('diversity_analyzer')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    params = parser.parse_args()

    CheckBinariesExistance(params, log)
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

    PrepareConfigs(params, log)
    try:
        cdr_command_line = run_cdr_labeler + " " + params.cdr_labeler_config_file
        support.sys_call(cdr_command_line, log)
        if not params.skip_plots:
            log.info("\n==== Visualization of diversity statistics ====")
            visualize_vj_stats.main(["", os.path.join(params.output_dir, "cdr_details.txt"),
                                 os.path.join(params.output_dir, "shm_details.txt"),
                                 params.output_dir, log])
            log.info("\n==== Annotation report creation ====")
            html_report_writer.main(os.path.join(params.output_dir, "cdr_details.txt"),
                                os.path.join(params.output_dir, "shm_details.txt"),
                                os.path.join(params.output_dir, "plots"),
                                os.path.join(params.output_dir, "annotation_report.html"), log)
        Cleanup(params, log)
        log.info("\nThank you for using " + tool_name + "!\n")
    except (KeyboardInterrupt):
         log.info("\n" + tool_name + " was interrupted!")
    except Exception:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            log.info("\nERROR: Exception caught.")
             #supportInfo(log)
    except BaseException:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            log.info("\nERROR: Exception caught.")
            #supportInfo(log)

    log.info("Log was written to " + params.log_filename)


if __name__ == '__main__':
    main(sys.argv)
