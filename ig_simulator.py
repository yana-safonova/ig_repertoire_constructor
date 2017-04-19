#!/usr/bin/env python2

import os
import sys
import init
import logging
import shutil
import ntpath

import process_cfg
import support

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
cdr_labeler_config_dir = os.path.join(home_directory, "configs", "cdr_labeler")
ig_simulator_config_dir = os.path.join(home_directory, "configs", "ig_simulator")
ig_simulator_bin = os.path.join(home_directory, "build", "release", "bin", "ig_simulator")
test_dir = "ig_simulator_test"

tool_name = "IgSimulator"


def CheckBinariesExistance(params, log):
    if not os.path.exists(ig_simulator_bin):
        log.info("ERROR: Binary files were not found. Please compile " + tool_name + " before running.")
        sys.exit(1)


def TreeStrategyCorrect(tree_strategy):
    return tree_strategy in ["uniform", "wide", "deep"]


def LociParamCorrect(loci):
    return loci in ["IGH", "IGK", "IGL"]


def CheckParamsCorrectness(params, log):
    if not LociParamCorrect(params.loci):
        log.info("Loci " + params.loci + " is not recognized")
        sys.exit(1)
    if not TreeStrategyCorrect(params.tree_strategy):
        log.info("Tree Strategy " + params.tree_strategy + " is not recognized")
        sys.exit(1)


def SetOutputParams(params, log):
    params.output_dir = os.path.abspath(params.output_dir)
    params.output_config_dir = os.path.join(params.output_dir, "configs")


def PrepareOutputDir(params):
    if os.path.exists(params.output_dir):
        shutil.rmtree(params.output_dir)
    os.mkdir(params.output_dir)


def parse_args():
    from argparse import ArgumentParser
    parser = ArgumentParser(description="== " + tool_name + ": a tool for simulating antibody repertoires, clonal lineages and trees ==",
                            epilog="In case you have troubles running " + tool_name + ", you can write to igtools_support@googlegroups.com."
                            "Please provide us with igsimulator.log file from the output directory.",
                            add_help=False)
    req_args = parser.add_argument_group("Required params")
    output_args = req_args.add_mutually_exclusive_group(required=True)
    output_args.add_argument("-o", "--output",
                             type=str,
                             default="",
                             dest="output_dir",
                             help="Output directory")

    output_args.add_argument("--test",
                             action="store_const",
                             const=test_dir,
                             dest="output_dir",
                             help="Running in test mode")

    optional_args = parser.add_argument_group("Optional arguments")

    optional_args.add_argument("-l", "--loci",
                               type=str,
                               default="IGH",
                               dest="loci",
                               help="Loci: IGH, IGK, IGL" # ", TRA, TRB, TRG, TRD, TR (all TCRs) or all. "
                                    "[default: %(default)s]")

    optional_args.add_argument("-s", "--tree_strategy",
                               type=str,
                               default="deep",
                               dest="tree_strategy",
                               help="Tree strategy to use: uniform, wide, deep [default: %(default)s]")

    optional_args.add_argument("-n", "--n_metaroots",
                               type=int,
                               default=1000,
                               dest="number_of_metaroots")

    optional_args.add_argument("-h", "--help",
                               action="help",
                               help="Help message and exit")
    return parser.parse_args()


def get_logger():
    log = logging.getLogger(tool_name)
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    return log


def add_log_handler(params, log):
    # log file
    params.log_filename = os.path.join(params.output_dir, "ig_simulator.log")
    if os.path.exists(params.log_filename):
        os.remove(params.log_filename)
    log_handler = logging.FileHandler(params.log_filename, mode='a')
    log.addHandler(log_handler)
    return log


def PrintParams(params, log):
    log.info(tool_name + " parameters:")
    log.info("  Output directory:\t" + params.output_dir + "\n")
    log.info("  Loci:\t\t\t" + params.loci)
    log.info("  # of metaroots:\t\t" + str(params.number_of_metaroots) + "\n")
    log.info("  Tree strategy:\t\t" + params.tree_strategy + "\n")

########################################################################################################################


def CopyConfigs(params, log):
    if os.path.exists(params.output_config_dir):
        shutil.rmtree(params.output_config_dir)
    params.cdr_labeler_config_dir = os.path.abspath(os.path.join(params.output_config_dir, "cdr_labeler"))
    shutil.copytree(ig_simulator_config_dir, params.output_config_dir)
    shutil.copytree(cdr_labeler_config_dir, params.cdr_labeler_config_dir)
    params.output_config_file = os.path.join(params.output_config_dir, "config.info")
    if not os.path.exists(params.output_config_file):
        log.info("ERROR: Config file " + params.output_config_file + " was not found")
        sys.exit(1)


def ModifyConfigFiles(params, log):
    igs_params_dict = dict()
    igs_params_dict['output_dir'] = params.output_dir
    igs_params_dict['loci'] = params.loci
    igs_params_dict['number_of_metaroots'] = params.number_of_metaroots
    igs_params_dict['pool_manager_strategy'] = params.tree_strategy

    process_cfg.substitute_params(params.output_config_file, igs_params_dict, log)


def PrepareConfigs(params, log):
    CopyConfigs(params, log)
    ModifyConfigFiles(params, log)


def RunTool(params, log):
    try:
        igs_command_line = ig_simulator_bin + " " + \
                           params.output_config_file
        support.sys_call(igs_command_line, log)
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
    log = get_logger()
    params = parse_args()
    print(params)
    CheckBinariesExistance(params, log)
    CheckParamsCorrectness(params, log)
    SetOutputParams(params, log)

    PrepareOutputDir(params)
    log = add_log_handler(params, log)

    # print command line
    command_line = "Command_line: "
    command_line += " ".join(argv)
    log.info(command_line + "\n")
    PrintParams(params, log)
    log.info("Log will be written to " + params.log_filename + "\n")

    PrepareConfigs(params, log)

    RunTool(params, log)
    log.info("Log was written to " + params.log_filename)

if __name__ == '__main__':
    main(sys.argv)
