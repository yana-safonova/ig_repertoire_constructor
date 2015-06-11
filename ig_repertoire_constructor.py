#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import os
import logging
import shutil
import argparse

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'
spades_src = os.path.join(home_directory, "src/python_pipeline/")
config_dir = os.path.join(home_directory, "configs/ig_repertoire_constructor/")
ig_binary = os.path.join(home_directory, "build/release/bin/ig_repertoire_constructor")

sys.path.append(spades_src)
import process_cfg

def supportInfo(log):
    log.info("\nIn case you have troubles running IgRepertoireConstructor, you can write to igtools_support@googlegroups.com.")
    log.info("Please provide us with ig_repertoire_constructor.log file from the output directory.")

def SetOutputParams(params, output_dir):
    params.output_dir = os.path.abspath(output_dir)
    params.dataset_file = os.path.join(params.output_dir, params.dataset_file)
    params.config_dir = os.path.join(params.output_dir, params.config_dir)
    params.config_file = os.path.join(params.config_dir, params.config_file)
    params.saves_dir = os.path.join(params.output_dir, params.saves_dir)
    params.temp_files_dir = os.path.join(params.output_dir, params.temp_files_dir)
    params.result_clusters = os.path.join(params.output_dir, "constructed_repertoire.clusters.fa")
    params.result_rcm = os.path.join(params.output_dir, "constructed_repertoire.rcm")

def PrepareOutputDir(params):
    if params.entry_point == "ig_repertoire_constructor" and os.path.exists(params.output_dir):
        shutil.rmtree(params.output_dir)
    if not os.path.isdir(params.output_dir):
        os.makedirs(params.output_dir)

def CheckParamsCorrectness(params, log, parser):
    if not "output" in params or params.output == "":
        log.info("ERROR: Output directory (-o) was not specified")
        parser.print_help()
        sys.exit(-1)
    if not "reads" in params or params.reads == "":
        log.info("ERROR: Reads (-s) were not specified")
        parser.print_help()
        sys.exit(-1)
    if not os.path.exists(params.reads):
        log.info("ERROR: File with reads " + params.reads + " were not found")
        parser.print_help()
        sys.exit(-1)
    if not os.path.isabs(params.reads):
        params.reads = os.path.abspath(params.reads)

def PrintParams(params, log):
    log.info("IgRepertoireConstructor parameters:")
    log.info("  Input reads:\t\t\t" + params.reads)
    log.info("  Output directory:\t\t" + params.output_dir)
    log.info("  Number of threads:\t\t" + str(params.num_threads))
    log.info("  Overlap mismatches threshold:\t" + str(params.mismatches_threshold))
    log.info("  Memory limit (in Gb):\t\t" + str(params.max_memory))
    log.info("  Entry point:\t\t\t" + params.entry_point)

def CreateDatasetYaml(params, log):
    dataset_fhandler = open(params.dataset_file, "w")
    dataset_fhandler.write("- single reads: [" + params.reads + "]\n")
    dataset_fhandler.write("  type: single\n")
    dataset_fhandler.close()
    log.info("Dataset was written to " + params.dataset_file)

def CopyConfigs(params, log):
    if os.path.exists(params.config_dir):
        shutil.rmtree(params.config_dir)
    shutil.copytree(config_dir, params.config_dir)

def CreateParamDict(params):
    param_dict = dict()
    param_dict['output_dir'] = params.output_dir
    param_dict['dataset'] = params.dataset_file
    param_dict['output_saves'] = params.saves_dir
    param_dict['temp_files'] = params.temp_files_dir
    param_dict['load_from'] = params.saves_dir
    param_dict['threads_count'] = params.num_threads
    param_dict['entry_point'] = params.entry_point
    param_dict['overlap_mismatches_threshold'] = params.mismatches_threshold
    param_dict['max_memory'] = params.max_memory
    param_dict['class_joining_edge_threshold'] = params.joint_thresh
    if params.output_dense_sgraphs:
        param_dict["output_dense_subgraphs"] = "true"
    return param_dict

def PrepareConfigs(params, log):
    CopyConfigs(params, log)
    param_dict = CreateParamDict(params)
    if not os.path.exists(params.config_file):
        log.info("ERROR: config file was not found")
        sys.exit(1)
    process_cfg.substitute_params(params.config_file, param_dict, log)

def CheckIgRepertoireConstructor(params, log):
    if os.path.exists(params.result_clusters):
        log.info("\n * CLUSTERS.FASTA for final repertoire is in " + params.result_clusters)
    else:
        log.info("\nERROR: CLUSTERS.FASTA for final repertoire was not found")
        supportInfo(log)
        raise SystemExit()
    if os.path.exists(params.result_rcm):
        log.info(" * RCM for final repertoire is in " + params.result_rcm)
    else:
        log.info("ERROR: RCM for final repertoire was not found")
        supportInfo(log)
        raise SystemExit()

def RunIgRepertoireConstructor(params, log):
    if not os.path.exists(ig_binary):
        log.info("\nERROR: IgRepertoireConstructor binary file was not found!")
        sys.exit(1)
    command_line = ig_binary + " " + params.config_file
    log.info("\n==== IgRepertoireConstructor starts\n")

    err_code = os.system(command_line + " 2>&1 | tee -a " + params.log_filename)
    if err_code != 0:
        log.info("\nERROR: IgRepertoireConstructor was finished abnormally, error code: " + str(err_code))
        supportInfo(log)
        sys.exit(-1)
    log.info("\n==== IgRepertoireConstructor finished\n")
    CheckIgRepertoireConstructor(params, log)
    log.info("\nThank you for using IgRepertoireConstructor!\n")

def CleanOutputDir(params, log):
    #log.info("Removing temporary data")
    shutil.rmtree(params.temp_files_dir)
    if not params.save_hamming_graphs:
        #log.info("Removing Hamming graphs")
        for fname in os.listdir(params.output_dir):
            path = os.path.join(params.output_dir, fname)
            if fname.startswith("hamming_graphs_tau_") and os.path.isdir(path):
                shutil.rmtree(path)

def main():
    # Parse commandline args
    parser = argparse.ArgumentParser(description="TODO Add some description",
                                     epilog="""
    In case you have troubles running IgRepertoireConstructor, you can write to igtools_support@googlegroups.com.
    Please provide us with ig_repertoire_constructor.log file from the output directory.
                                     """,
                                     add_help=False)

    req_args = parser.add_argument_group("Input")
    input_args = req_args.add_mutually_exclusive_group(required=True)
    input_args.add_argument("-s", "--reads",
                            type=str,
                            default="", # FIXME This is only for ace's version of python. Locally it works great w/o it
                            help="cleaned FASTQ reads corresponding to variable regions of immunoglobulins")
    input_args.add_argument("--test",
                            action="store_const",
                            const="test_dataset/merged_reads.fastq",
                            dest="reads",
                            help="`merged_reads` test dataset")
    input_args.add_argument("--testIGHV",
                            action="store_const",
                            const="test_dataset/IGHV1-8.fastq",
                            dest="reads",
                            help="IGHV1-8 test dataset")
    input_args.add_argument("--test2",
                            action="store_const",
                            const="test_dataset/test2.fastq",
                            dest="reads",
                            help="test dataset based on 2_SAM13306970 data using 13 largest connectivity components")
    input_args.add_argument("--test7",
                            action="store_const",
                            const="test_dataset/test7.fastq",
                            dest="reads",
                            help="test dataset based on 7_SAM15574987 data using 13 largest connectivity components. Be careful, it's longtime")

    out_args = parser.add_argument_group("Output")
    out_args.add_argument("-o", "--output",
                          type=str,
                          default="ig_repertoire_constructor_test",
                          help="output directory [default \"%(default)s\"]")

    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument("-t", "--threads",
                               type=int,
                               default=16,
                               dest="num_threads",
                               help="threads number [default %(default)d]")
    optional_args.add_argument("-m", "--memory",
                               type=int,
                               default=250,
                               dest="max_memory",
                               help="RAM limit for ig_repertoire_constructor in Gb [default: %(default)d]")
    optional_args.add_argument("--tau",
                               type=int,
                               default=3,
                               dest="mismatches_threshold",
                               help="maximum allowed mismatches between reads in cluster [default: %(default)d]")
    optional_args.add_argument("-h", "--help",
                               action="help",
                               help="show this help message and exit")

    dev_args = parser.add_argument_group("Developer arguments")
    dev_args.add_argument("--joint-thresh",
                          type=float,
                          default=0.6,
                          help="threshold for minimum value of edge fill-in in dense subgraph construction procedure [default %(default)2.1f]")
    dev_args.add_argument('--entry-point',
                          type=str,
                          default="ig_repertoire_constructor",
                          help="continue from the given stage [default %(default)s]")
    shg_args = dev_args.add_mutually_exclusive_group(required=False)
    shg_args.add_argument("--save-hgraphs",
                          action="store_const",
                          const=True,
                          dest="save_hamming_graphs",
                          help="saves Hamming graphs in GRAPH format [default]")
    shg_args.add_argument("--no-save-hgraphs",
                          action="store_const",
                          const=False,
                          dest="save_hamming_graphs",
                          help="")
    ods_args = dev_args.add_mutually_exclusive_group(required=False)
    ods_args.add_argument("--output-dense-sgraphs",
                          action="store_const",
                          const=True,
                          dest="output_dense_sgraphs",
                          help="outputs decomposition into dense subgraphs [default]")
    ods_args.add_argument("--no-output-dense-sgraphs",
                          action="store_const",
                          const=False,
                          dest="output_dense_sgraphs",
                          help="")

    parser.set_defaults(output_dense_sgraphs=True,
                        save_hamming_graphs=True)

    parser.set_defaults(dataset_file="dataset.yaml",
                        config_dir="configs",
                        config_file="config.info",
                        saves_dir="saves",
                        temp_files_dir="temp_files")

    # prepare log
    log = logging.getLogger('ig_repertoire_constructor')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    # parse command line
    params = parser.parse_args()

    # params check
    CheckParamsCorrectness(params, log, parser)

    SetOutputParams(params, params.output)
    PrepareOutputDir(params)

    # Param print
    PrintParams(params, log)

    # log file
    params.log_filename = os.path.join(params.output_dir, "ig_repertoire_constructor.log")
    if os.path.exists(params.log_filename):
        os.remove(params.log_filename)
    log_handler = logging.FileHandler(params.log_filename, mode='a')
    log.addHandler(log_handler)
    log.info("Log will be written to " + params.log_filename + "\n")

    # print command line
    command_line = "Command line: " + " ".join(sys.argv)
    log.info(command_line)

    # prepare dataset.yaml
    CreateDatasetYaml(params, log)

    # prepare configs
    PrepareConfigs(params, log)

    try:
        # run IgRepertoireConstructor
        RunIgRepertoireConstructor(params, log)
    except (KeyboardInterrupt):
        log.info("\nIgRepertoireConstructor was interrupted!")
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

    CleanOutputDir(params, log)

    log.info("Log was written to " + params.log_filename)

if __name__ == '__main__':
    main()
