#!/usr/bin/env python2

import os
import sys
import shutil
import logging


home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
shm_kmer_matrix_estimator_config_dir = \
        os.path.join(home_directory, "configs", "shm_kmer_matrix_estimator")
shm_kmer_matrix_estimator_bin = \
        os.path.join(home_directory, "build", "release", "bin",
                     "shm_kmer_matrix_estimator")

py_src = os.path.join(home_directory, "src/python_pipeline/")
sys.path.append(py_src)
import process_cfg
import support

tool_name = "SHM Kmer Matrix Estimator"


def CheckBinariesExistance(log):
    if not os.path.exists(shm_kmer_matrix_estimator_bin):
        log.info("ERROR: Binary files were not found. Please compile " +
                 tool_name + " before running.")
        sys.exit(1)


def PrepareParams(params, log):
    params.v_alignments = os.path.abspath(params.v_alignments)
    params.cdr_details = os.path.abspath(params.cdr_details)
    params.output_dir = os.path.abspath(params.output_dir)
    params.output_config_dir = os.path.join(params.output_dir, "configs")
    params.output_config_dir = os.path.abspath(params.output_config_dir)


def PrepareOutputDir(params, log):
    if not os.path.exists(params.output_dir):
        # log.info("Cleaning output dir")
        # shutil.rmtree(params.output_dir)
        os.mkdir(params.output_dir)


def PrintParams(params, log):
    log.info(tool_name + " parameters:")
    log.info("  Input V alignemnts:\t\t" + params.v_alignments)
    log.info("  Input CDR detailss:\t\t" + params.cdr_details)
    log.info("  Output directory:\t\t" + params.output_dir + "\n")


def PrepareLog():
    log = logging.getLogger('diversity_analyzer')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    return log


def PrepareParser():
    def extant_file(x):
        import argparse
        """
        'Type' for argparse - checks that file exists but does not open.
        """
        if not os.path.exists(x):
            raise argparse.ArgumentTypeError("%s does not exist" % x)
        return x

    from src.python_add.argparse_ext import ArgumentHiddenParser
    parser = ArgumentHiddenParser(description="== " + tool_name + " ==",
                                  add_help=False)
    req_args = parser.add_argument_group("Required params")
    req_args.add_argument("-v", "--v_alignments",
                          dest="v_alignments",
                          type=extant_file,
                          required=True,
                          help="Input V alignments FASTA")
    req_args.add_argument("-c", "--cdr_details",
                          type=extant_file,
                          dest="cdr_details",
                          required=True,
                          help="Input CDR Details")
    req_args.add_argument("-o", "--output",
                          type=str,
                          dest="output_dir",
                          required=True,
                          help="Output directory")

    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument("-s", "--mutation_strategy",
                               type=str,
                               default="NoKNeighbours",
                               dest="mutation_strategy",
                               choices=["Trivial", "NoKNeighbours"],
                               help="Mutation Strategy for considering " +
                                    "close mutations: [default: %(default)s]")
    optional_args.add_argument("-f", "--functionality",
                               type=str,
                               default="all",
                               dest="functionality_method",
                               choices=["all", "productive", "nonproductive"],
                               help="Functionality Strategy : [default: %(default)]")
    optional_args.add_argument("-h", "--help",
                               action="help",
                               help="Help message and exit")
    return parser


def PrepareLogFilename(params, log):
    params.log_filename = os.path.join(params.output_dir,
                                       "shm_kmer_matrix_estimator.log")
    if os.path.exists(params.log_filename):
        os.remove(params.log_filename)
    log_handler = logging.FileHandler(params.log_filename, mode='a')
    log.addHandler(log_handler)


def PrintCommandLine(argv, params, log):
    command_line = "Command_line: "
    command_line += " ".join(argv)
    log.info(command_line + "\n")
    PrintParams(params, log)
    log.info("Log will be written to " + params.log_filename + "\n")


###############################################################################
def CopyConfig(params, log):
    if os.path.exists(params.output_config_dir):
        shutil.rmtree(params.output_config_dir)
    shutil.copytree(shm_kmer_matrix_estimator_config_dir,
                    params.output_config_dir)
    params.output_config_file = \
        os.path.join(params.output_config_dir, "configs.info")


def ModifyConfigFiles(params, log):
    shm_param_dict = dict()
    shm_param_dict["v_alignments"] = params.v_alignments
    shm_param_dict["cdr_details"] = params.cdr_details
    shm_param_dict["output_dir"] = params.output_dir
    shm_param_dict["mutations_strategy_method"] = params.mutation_strategy
    shm_param_dict["functionality_method"] = params.functionality_method
    process_cfg.substitute_params(params.output_config_file,
                                  shm_param_dict, log)


def PrepareConfig(params, log):
    CopyConfig(params, log)
    ModifyConfigFiles(params, log)

###############################################################################


def RunTool(params, log):
    try:
        cdr_command_line = shm_kmer_matrix_estimator_bin + " " + \
                           params.output_config_file
        support.sys_call(cdr_command_line, log)
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
    except BaseException:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            log.info("\nERROR: Exception caught.")


def main(argv):
    log = PrepareLog()
    CheckBinariesExistance(log)
    parser = PrepareParser()
    params = parser.parse_args()

    PrepareParams(params, log)
    PrepareOutputDir(params, log)
    PrepareLogFilename(params, log)
    PrintCommandLine(argv, params, log)

    PrepareConfig(params, log)

    RunTool(params, log)
    log.info("Log was written to " + params.log_filename)


if __name__ == '__main__':
    main(sys.argv)
