#!/usr/bin/env python2

import sys
import os
import init
import logging
import shutil
import ntpath

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'
pipeline_dir = os.path.join(home_directory, "py/pipeline/")
config_dir = os.path.join(home_directory, "configs/dense_sgraph_finder/")
ig_binary = os.path.join(home_directory, "build/release/bin/ig_repertoire_constructor")

sys.path.append(pipeline_dir)
import process_cfg
import support

def CheckParamsCorrectness(params, log, parser):
    if not "output" in params or params.output == "":
        log.info("ERROR: Output directory (-o) was not specified")
        parser.print_help()
        sys.exit(-1)
    if not "graph" in params or params.graph == "":
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
    log.info("  Input graph:\t\t\t\t\t" + params.graph)
    log.info("  Output directory:\t\t\t\t" + params.output)
    log.info("  Number of threads:\t\t\t\t" + str(params.num_threads))
    log.info("  Minimum size of processed graphs:\t\t" + str(params.min_graph_size))
    log.info("  Minimum edge fill-in of dense subgraph:\t" + str(params.min_fillin))
    if params.create_trivial_decomposition:
        log.info("  Trivial decomposition will be created")
    log.info("  Output auxiliary files:\t\t\t" + str(params.save_aux_files))
    log.info("\n")

def supportInfo(log):
    log.info("\nIn case you have troubles running DSF, you can write to igtools_support@googlegroups.com.")
    log.info("Please provide us with dense_sungraph_finder.log file from the output directory.")

def SetOutputParams(params, output_dir):
    params.output = os.path.abspath(output_dir)
    params.config_dir = os.path.join(params.output, params.config_dir)
    params.config_file = os.path.join(params.config_dir, params.config_file)

def PrepareOutputDir(params, log):
    log.info("Preparing output dir")
    Cleanup(params, log)
    if params.clean_output_dir and os.path.exists(params.output):
        log.info("Removing %s" % params.output)
        shutil.rmtree(params.output)
    if not os.path.isdir(params.output):
        os.makedirs(params.output)

def CopyConfigs(params, log):
    if os.path.exists(params.config_dir):
        log.info("Removing %s" % params.config_dir)
        shutil.rmtree(params.config_dir)
    shutil.copytree(config_dir, params.config_dir)

def CreateParamDict(params):
    param_dict = dict()
    param_dict['output_dir'] = params.output
    param_dict['graph_filename'] = params.graph
    param_dict['threads_count'] = params.num_threads
    param_dict['min_fillin_threshold'] = params.min_fillin
    param_dict['min_graph_size'] = params.min_graph_size
    param_dict['create_trivial_decomposition'] = process_cfg.bool_to_str(params.create_trivial_decomposition)
    param_dict['path_to_metis'] = os.path.join(home_directory, "build/release/bin/")
    param_dict['min_supernode_size'] = params.min_snode_size
    return param_dict

def PrepareConfigs(params, log):
    CopyConfigs(params, log)
    param_dict = CreateParamDict(params)
    if not os.path.exists(params.config_file):
        log.info("ERROR: config file was not found")
        sys.exit(1)
    process_cfg.substitute_params(params.config_file, param_dict, log)

def Cleanup(params, log):
    log.info("Cleaning up")
    params.sgraph_dir = os.path.join(params.output, "connected_components")
    params.decomposition_dir = os.path.join(params.output, "dense_subgraphs")
    params.metis_output = os.path.join(params.output, "metis.output")
    if os.path.exists(params.metis_output):
        log.info("Removing %s" % params.metis_output)
        os.remove(params.metis_output)
    if not "save_aux_files" in params or not params.save_aux_files:
        if os.path.exists(params.sgraph_dir):
            log.info("Removing %s" % params.sgraph_dir)
            shutil.rmtree(params.sgraph_dir)
        if os.path.exists(params.decomposition_dir):
            log.info("Removing %s" % params.decomposition_dir)
            shutil.rmtree(params.decomposition_dir)

def main(argv, external_logger = ""):
    from argparse import ArgumentParser
    parser = ArgumentParser(description="== DSF: an algorithm for corrupted cliques search ==",
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
                            dest="graph",
                            help="Input graph in GRAPH format")
    input_args.add_argument("--test",
                            action="store_const",
                            const=os.path.join(home_directory, "test_dataset/dsf/test.graph"),
                            dest="graph",
                            help="Running test dataset")

    out_args = parser.add_argument_group("Output")
    out_args.add_argument("-o", "--output",
                            type=str,
                            default=os.path.join(home_directory, "dsf_test"),
                            help="Output directory")

    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument("-t", "--threads",
                               type=int,
                               default=16,
                               dest="num_threads",
                               help="Threads number [default: %(default)d]")
    optional_args.add_argument("-f", '--min-fillin',
                               type=float,
                               default=0.6,
                               dest="min_fillin",
                               help='Minimum fill-in of dense subgraphs [default: %(default)f]')
    optional_args.add_argument("-n", "--min-snode-size",
                               type=int,
                               default=5,
                               dest="min_snode_size",
                               help="Minimum vertex weight that prevents its gluing with other heavy vertex "
                                    "[default: %(default)d]")
    optional_args.add_argument("-s", "--min-size",
                               type=int,
                               default=5,
                               dest="min_graph_size",
                               help="Minimum size of graph where dense subgraphs will be computed "
                                    "[default: %(default)d]")
    optional_args.add_argument("--create-triv-dec",
                               action="store_const",
                               const=True,
                               dest="create_trivial_decomposition",
                               help='Creating decomposition according to connected components [default: False]')
    optional_args.add_argument("--save-aux-files",
                               action="store_const",
                               const=True,
                               dest="save_aux_files",
                               help="Saving auxiliary files: subgraphs in GRAPH format and their decompositions "
                                    "[default: False]")
    optional_args.add_argument("--clean-output-dir",
                               default=True,
                               dest="clean_output_dir",
                               action="store_true",
                               help="Clean output directory on start [default]")
    optional_args.add_argument("--no-clean-output-dir",
                               default=True,
                               dest="clean_output_dir",
                               action="store_false",
                               help="Do not clean output directory on start")
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
    if external_logger != "":
        external_log_handler = logging.FileHandler(external_logger, mode = "a")
        log.addHandler(external_log_handler)

    args = [arg for arg in argv if ntpath.basename(arg) != 'dense_subgraph_finder.py']
    params = parser.parse_args(args)

    CheckParamsCorrectness(params, log, parser)
    SetOutputParams(params, params.output)

    PrepareOutputDir(params, log)

    # log file
    params.log_filename = os.path.join(params.output, "dense_subgraph_finder.log")
    if os.path.exists(params.log_filename):
        log.info("Removing %s" % params.log_filename)
        os.remove(params.log_filename)
    log_handler = logging.FileHandler(params.log_filename, mode='a')
    log.addHandler(log_handler)

    # print command line
    command_line = "Command_line: "
    if argv[0] != "dense_subgraph_finder.py":
        command_line += "dense_subgraph_finder.py "
    command_line += " ".join(argv)
    log.info(command_line + "\n")
    PrintParams(params, log)
    log.info("Log will be written to " + params.log_filename + "\n")

    PrepareConfigs(params, log)

    # run dense subgraph finder
    try:
        dsf_command_line = init.PathToBins.run_dense_sgraph_finder + " " + params.config_file
        support.sys_call(dsf_command_line, log)
        Cleanup(params, log)
        log.info("\nThank you for using Dense Subgraph Finder!\n")
    except (KeyboardInterrupt):
        log.info("\nDense subgraph finder was interrupted!")
    except Exception:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            log.info("\nERROR: Exception caught.")
            supportInfo(log)
    except BaseException:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            log.info("\nERROR: Exception caught.")
            supportInfo(log)

    log.info("Log was written to " + params.log_filename)


if __name__ == '__main__':
    main(sys.argv)
