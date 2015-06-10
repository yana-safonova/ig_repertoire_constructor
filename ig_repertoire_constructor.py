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

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'
spades_src = os.path.join(home_directory, "src/python_pipeline/")
config_dir = os.path.join(home_directory, "configs/ig_repertoire_constructor/")
ig_binary = os.path.join(home_directory, "build/release/bin/ig_repertoire_constructor")

sys.path.append(spades_src)
import process_cfg

def usage(log, show_hidden=False):
    log.info("./ig_repertoire_constructor.py [options] -s <filename> -o <output_dir>")
    log.info("\nBasic options:")
    log.info("  -s\t\t\t<filename>\tcleaned FASTQ reads corresponding to variable regions of immunoglobulins (required)")
    log.info("  -o/--output\t\t<output_dir>\toutput directory (required)")
    log.info("  -t/--threads\t\t<int>\t\tthreads number [default: 16]")
    log.info("  --test\t\t\t\truns test dataset")
    log.info("  --help\t\t\t\tprints help")

    log.info("\nAdvanced options:")
    log.info("  --tau\t\t\t<int>\t\tmaximum allowed mismatches between reads in cluster [default: 3]")
    log.info("  -m/--memory\t\t<int>\t\tRAM limit for ig_repertoire_constructor in Gb [default: 250]")

    if show_hidden:
        log.info("\nHidden options:")
        log.info("  --entry-point\t\t<stage_name>\tcontinue from the given stage")
        log.info("  --help-hidden\t\t\t\tprints this usage message with all hidden options")
        log.info("  --save-hgraphs\t\t\tsaves Hamming graphs in GRAPH format")
        log.info("  --output-dense-sgraphs\t\toutputs decomposition into dense subgraphs")

def supportInfo(log):
    log.info("In case you have troubles running IgRepertoireConstructor, you can write to igtools_support@googlegroups.com.")
    log.info("Please provide us with ig_repertoire_constructor.log file from the output directory.")

def SetOutputParams(params, output_dir):
    params.output_dir = os.path.abspath(output_dir)
    params.dataset_file = os.path.join(params.output_dir, params.dataset_file)
    params.config_dir = os.path.join(params.output_dir, params.config_dir)
    params.config_file = os.path.join(params.config_dir, params.config_file)
    params.saves_dir = os.path.join(params.output_dir, params.saves_dir)
    params.temp_files_dir = os.path.join(params.output_dir, params.temp_files_dir)

def PrepareOutputDir(output_dir):
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

def CheckParamsCorrectness(params, log):
    if params.output_dir == "":
        log.info("ERROR: Output directory (-o) was not specified")
        usage(log)
        sys.exit(-1)
    if params.reads == "":
        log.info("ERROR: Reads (-s) were not specified")
        usage(log)
        sys.exit(-1)
    if not os.path.exists(params.reads):
        log.info("ERROR: File with reads " + params.reads + " were not found")
        usage(log)
        sys.exit(-1)

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
    log.info("\n * CLUSTERS.FASTA for final repertoire is in " + params.output_dir + "/constructed_repertoire.clusters.fa")
    log.info(" * RCM for final repertoire is in " + params.output_dir + "/constructed_repertoire.rcm")
    log.info("\nThank you for using IgRepertoireConstructor!\n")

def CleanOutputDir(params, log):
    log.info("Removing temporary data")
    shutil.rmtree(params.temp_files_dir)
    if not params.save_hamming_graphs:
        log.info("Removing Hamming graphs")
        for fname in os.listdir(params.output_dir):
            path = os.path.join(params.output_dir, fname)
            if fname.startswith("hamming_graphs_tau_") and os.path.isdir(path):
                shutil.rmtree(path)

def main():
    import argparse
    # Parse commandline args
    parser = argparse.ArgumentParser(description="TODO Add some description",
                                     epilog="""
    In case you have troubles running IgRepertoireConstructor, you can write to igtools_support@googlegroups.com.
    Please provide us with ig_repertoire_constructor.log file from the output directory.
                                     """)

    parser.add_argument("-o", "--output",
                        type=str,
                        # required=True,
                        help="output directory (required)")
    parser.add_argument("-s", "--reads",
                        type=str,
                        # required=True,
                        help="cleaned FASTQ reads corresponding to variable regions of immunoglobulins (required)")
    parser.add_argument("-t", "--threads",
                        type=int,
                        default=16,
                        dest="num_threads",
                        help="threads number [default %(default)d]")
    parser.add_argument('--entry-point',
                        type=str,
                        default="ig_repertoire_constructor",
                        help="continue from the given stage")
    parser.add_argument("--tau",
                        type=int,
                        default=3,
                        dest="mismatches_threshold",
                        help="maximum allowed mismatches between reads in cluster [default: %(default)d]")
    parser.add_argument("-m", "--memory",
                        type=int,
                        default=250,
                        dest="max_memory",
                        help="RAM limit for ig_repertoire_constructor in Gb [default: %(default)d]")
    parser.add_argument("--joint-thresh",
                        type=float,
                        default=0.3,
                        help="Yana, please, provide some description;)")
    parser.add_argument("--save-hgraphs",
                        action="store_true",
                        dest="save_hamming_graphs",
                        help="outputs decomposition into dense subgraphs")
    parser.add_argument("--output-dense-sgraphs",
                        action="store_true",
                        help="saves Hamming graphs in GRAPH format")
    parser.add_argument("--test",
                        action="store_true",
                        help="runs test dataset")

    # prepare log
    log = logging.getLogger('ig_repertoire_constructor')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    # parse command line
    params = parser.parse_args()
    params.dataset_file = "dataset.yaml"
    params.config_dir = "configs"
    params.config_file = "config.info"
    params.saves_dir = "saves"
    params.temp_files_dir = "temp_files"
    if params.test:
        SetOutputParams(params, os.path.abspath('ig_repertoire_constructor_test'))
        params.reads = os.path.join(home_directory, 'test_dataset/merged_reads.fastq')
        params.reads = os.path.abspath(params.reads)
    else:
        SetOutputParams(params, params.output)
    if not os.path.isabs(params.reads):
        params.reads = os.path.abspath(params.reads)

    PrepareOutputDir(params.output_dir)

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

    # params check and print
    CheckParamsCorrectness(params, log)
    PrintParams(params, log)

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
