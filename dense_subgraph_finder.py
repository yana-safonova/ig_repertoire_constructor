#!/usr/bin/env python
from cvxopt.modeling import op

import sys
import os
import init
import logging
import shutil

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'
spades_src = os.path.join(home_directory, "src/python_pipeline/")
config_dir = os.path.join(home_directory, "configs/dense_sgraph_finder/")
ig_binary = os.path.join(home_directory, "build/release/bin/ig_repertoire_constructor")

sys.path.append(spades_src)
import process_cfg

def CheckParamsCorrectness(params, log, parser):
    if not "output" in params or params.output == "":
        log.info("ERROR: Output directory (-o) was not specified")
        parser.print_help()
        sys.exit(-1)
    if not "graph" in params or params.reads == "":
        log.info("ERROR: Graph (-g/--graph) was not specified")
        parser.print_help()
        sys.exit(-1)
    if not os.path.exists(params.graph):
        log.info("ERROR: File with graph " + params.graph + " was not found")
        parser.print_help()
        sys.exit(-1)
    if not os.path.isabs(params.graph):
        params.graph = os.path.abspath(params.graph)

def PrintParams(params, log):
    log.info("DSF parameters:")
    log.info("  Input graph:\t\t\t" + params.graph)
    log.info("  Output directory:\t\t" + params.output)
    log.info("  Number of threads:\t\t" + str(params.num_threads))

def supportInfo(log):
    log.info("\nIn case you have troubles running DSF, you can write to igtools_support@googlegroups.com.")
    log.info("Please provide us with dense_sungraph_finder.log file from the output directory.")

def SetOutputParams(params, output_dir):
    params.output = os.path.abspath(output_dir)
    params.config_dir = os.path.join(params.output, params.config_dir)
    params.config_file = os.path.join(params.config_dir, params.config_file)

def PrepareOutputDir(params, log):
    if os.path.exists(params.output):
        shutil.rmtree(params.output)
    if not os.path.isdir(params.output):
        os.makedirs(params.output)

def CopyConfigs(params, log):
    if os.path.exists(params.config_dir):
        shutil.rmtree(params.config_dir)
    shutil.copytree(config_dir, params.config_dir)

def CreateParamDict(params):
    param_dict = dict()
    param_dict['output_dir'] = params.output
    param_dict['graph_filename'] = params.graph
    param_dict['threads_count'] = params.num_threads
    return param_dict

def PrepareConfigs(params, log):
    CopyConfigs(params, log)
    param_dict = CreateParamDict(params)
    if not os.path.exists(params.config_file):
        log.info("ERROR: config file was not found")
        sys.exit(1)
    process_cfg.substitute_params(params.config_file, param_dict, log)

def main():
    from src.python_add.argparse_ext import ArgumentHiddenParser
    parser = ArgumentHiddenParser(description="== DSF: an algorithm for corrupted cliques search ==",
                                  epilog="""
                                  In case you have troubles running DSF, you can write to igtools_support@googlegroups.com.
                                  Please provide us with dense_subgraph_finder.log file from the output directory.
                                  """,
                                  add_help=False)
    req_args = parser.add_argument_group("Input")
    input_args = req_args.add_mutually_exclusive_group(required=True)
    input_args.add_argument("-g", "--graph",
                            type=str,
                            default="",
                            help="Input graph in GRAPH format")

    out_args = parser.add_argument_group("Output")
    out_args.add_argument("-o", "--output",
                            type=str,
                            default="",
                            help="Output directory")

    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument("-t", "--threads",
                               type=int,
                               default=16,
                               dest="num_threads",
                               help="Threads number [default %(default)d]")
    optional_args.add_argument("--test",
                            action="store_const",
                            const="test_dataset/test.graph",
                            dest="reads",
                            help="Running test dataset")
    optional_args.add_argument("-h", "--help",
                               action="help",
                               help="Help message and exit")

    parser.set_defaults(config_dir="configs",
                        config_file="config.info")

    # prepare log
    log = logging.getLogger('dense_subgraph_finder')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    # parse command line
    params = parser.parse_args()

    CheckParamsCorrectness(params, log, parser)
    SetOutputParams(params, params.output)
    PrintParams(params, log)

    PrepareOutputDir(params, log)

    # log file
    params.log_filename = os.path.join(params.output, "dense_subgraph_finder.log")
    if os.path.exists(params.log_filename):
        os.remove(params.log_filename)
    log_handler = logging.FileHandler(params.log_filename, mode='a')
    log.addHandler(log_handler)
    log.info("Log will be written to " + params.log_filename + "\n")

    # print command line
    command_line = "Command line: " + " ".join(sys.argv)
    log.info(command_line)

    PrepareConfigs(params, log)

    # just draft of main fuction
    dsf_command_line = init.PathToBins.run_dense_sgraph_finder + " " + params.config_file
    print "Command line: " + dsf_command_line
    error_code = os.system(dsf_command_line)
    if error_code != 0:
        print "ERROR: Dense sgraph finder finished abnormally"
        sys.exit(1)

if __name__ == '__main__':
    main() 
