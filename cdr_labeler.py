#!/usr/bin/env python2

import os
import sys
import init
import logging
import shutil
import ntpath

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'
py_src = os.path.join(home_directory, "src/python_pipeline/")
cdr_labeler_config_dir = os.path.join(home_directory, "configs/cdr_labeler/")
vj_finder_config_dir = os.path.join(home_directory, "configs/vj_finder")
cdr_labeler_bin = "build/release/bin/cdr_labeler"
run_cdr_labeler = "build/release/bin/./cdr_labeler"

sys.path.append(py_src)
import process_cfg
import support

test_reads = os.path.join(home_directory, "test_dataset/merged_reads.fastq")
test_dir = os.path.join(home_directory, "cdr_test")

def CheckBinariesExistance(params, log):
    if not os.path.exists(cdr_labeler_bin):
        log.info("ERROR: CDR labeler binary file was not found. Please compile CDR labeler before running.")
        sys.exit(1)

def DomainParamCorrect(domain_str):
    return domain_str == "imgt" or domain_str == "kabat"

def LociParamCorrect(loci_str):
    return loci_str == "IG" or loci_str == "IGH" or loci_str == "IGK" or loci_str == "IGL"

def CheckParamsCorrectness(params, log):
    if not os.path.exists(params.input_reads):
        log.info("Input reads " + params.input_reads + " were not found")
        sys.exit(1)
    if not DomainParamCorrect(params.domain_system):
        log.info("Domain system " + params.domain_system + " is not recognized")
        sys.exit(1)
    if not LociParamCorrect(params.loci):
        log.info("Loci " + params.loci + " is not recognized")
        sys.exit(1)

def SetOutputParams(params, log):
    if params.input_reads == test_reads:
        params.output_dir = test_dir
    if params.input_reads != test_reads and params.output_dir == "":
        log.info("ERROR: Output dir (-o) was not specified")
        sys.exit(1)
    params.config_dir = os.path.join(params.output_dir, "configs")
    params.cdr_config_file = os.path.join(cdr_labeler_config_dir, "config.info")
    params.vj_finder_config_file = os.path.join(vj_finder_config_dir, "config.info")

def PrepareOutputDir(params):
    if os.path.exists(params.output_dir):
        shutil.rmtree(params.output_dir)
    os.mkdir(params.output_dir)

def PrintParams(params, log):
    log.info("CDR Labeler parameters:")
    log.info("  Input reads:\t\t" + params.input_reads)
    log.info("  Output directory:\t" + params.output_dir + "\n")
    log.info("  Loci:\t\t\t" + params.loci)
    log.info("  Domain system:\t" + params.domain_system + "\n")

########################################################################################################################

def CopyConfigs(params, log):
    if os.path.exists(params.config_dir):
        shutil.rmtree(params.config_dir)
    os.mkdir(params.config_dir)
    shutil.copytree(cdr_labeler_config_dir, os.path.join(params.config_dir, "cdr_labeler"))
    shutil.copytree(vj_finder_config_dir, os.path.join(params.config_dir, "vj_finder"))
    params.vj_finder_config_file = os.path.join(params.config_dir, "vj_finder/config.info")
    if not os.path.exists(params.vj_finder_config_file):
        log.info("ERROR: Config file " + params.vj_finder_config_file + " was not found")
        sys.exit(1)
    params.cdr_labeler_config_file = os.path.join(params.config_dir, "cdr_labeler/config.info")
    if not os.path.exists(params.cdr_labeler_config_file):
        log.info("ERROR: Config file " + params.cdr_labeler_config_file + " was not found")
        sys.exit(1)

def ModifyConfigFiles(params, log):
    cdr_param_dict = dict()
    cdr_param_dict['input_reads'] = params.input_reads
    cdr_param_dict['output_dir'] = params.output_dir
    cdr_param_dict['vj_finder_config'] = params.vj_finder_config_file
    cdr_param_dict['num_threads'] = params.num_threads
    cdr_param_dict['domain_system'] = params.domain_system
    process_cfg.substitute_params(params.cdr_labeler_config_file, cdr_param_dict, log)

    vj_param_dict = dict()
    vj_param_dict['loci'] = params.loci
    process_cfg.substitute_params(params.vj_finder_config_file, vj_param_dict, log)

def PrepareConfigs(params, log):
    CopyConfigs(params, log)
    ModifyConfigFiles(params, log)

########################################################################################################################

def main(argv):
    from src.python_add.argparse_ext import ArgumentHiddenParser
    parser = ArgumentHiddenParser(description="== CDR Labeler: a CDR search in full-length immunosequencing reads ==",
                                  epilog="""
                                  In case you have troubles running CDR Labeler, you can write to igtools_support@googlegroups.com.
                                  Please provide us with cdr_labeler.log file from the output directory.
                                  """,
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
                               help="Threads number [default: %(default)d]")
    optional_args.add_argument("-d", '--domain',
                               type=str,
                               default="imgt",
                               dest="domain_system",
                               help='Domain system for CDR search: imgt OR kabat [default: %(default)s]')

    vj_finder_args= parser.add_argument_group("VJ alignment params")
    optional_args.add_argument("-l", "--loci",
                               type=str,
                               default="IG",
                               dest="loci",
                               help="Loci: IGH, IGK, IGL, IG (all BCRs)" # ", TRA, TRB, TRG, TRD, TR (all TCRs) or all. "
                                    "[default: %(default)s]")

    optional_args.add_argument("-h", "--help",
                               action="help",
                               help="Help message and exit")

    #parser.set_defaults(config_dir=os.path.join(params.output_dir, "configs"),
    #                    cdr_config_file=os.path.join(cdr_labeler_config_dir, "config.info"),
    #                    vj_finder_config_file = os.path.join(vj_finder_config_dir, "config.info")))

    # prepare log
    log = logging.getLogger('cdr_labeler')
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
    params.log_filename = os.path.join(params.output_dir, "cdr_labeler.log")
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
        #Cleanup(params, log)
        log.info("\nThank you for using CDR Labeler!\n")
    except (KeyboardInterrupt):
         log.info("\nCDR Labeler was interrupted!")
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