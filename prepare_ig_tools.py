#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import getopt
import os
import logging
import shutil
import subprocess

import ig_tools_init

sys.path.append(ig_tools_init.spades_py_scripts_directory)

import process_cfg

class BaseOptions:
    long_options = "path-to-igblast= help dev-mode".split()
    short_options = ""

class Params:
    path_to_igblast = os.path.join(os.path.abspath(os.path.dirname(__file__)), \
        "src/tools/igblast/")

def usage(log):
    log.info("./prepare_ig_tools.py --path-to-igblast <path>")
    log.info("  --path-to-igblast\t<path>\t\t\tpath to IgBLAST tool [default: " + Params.path_to_igblast + "]")

def PrintParams(params, log):
    log.info("Path to IgBLAST tool " + params.path_to_igblast)

def CheckForIgblastCorrectness(path_to_blast, log):
    if os.path.exists(path_to_blast):
        log.info("* IgBlast tool will be run from " + os.path.abspath(path_to_blast))
        return
    log.info("ERROR: IgBlast was not found. Please use option --path-to-blast to specify it")
    usage(log)
    sys.exit(1)

def CheckForParamsCorrectness(params, log):
    log.info("\n======== Checking for IgBLAST paths correctness starts\n")
    CheckForIgblastCorrectness(params.path_to_igblast, log)
    log.info("\n======== Checking for IgBLAST paths correctness ends\n")

def PrepareConfig(params, path_to_config, log):
    params_dict = dict()
    params_dict["path_to_igblast"] = params.path_to_igblast
    process_cfg.substitute_params(path_to_config, params_dict, log)

def PrepareBinDirectory():
    if os.path.exists(ig_tools_init.ig_bin_directory):
        shutil.rmtree(ig_tools_init.ig_bin_directory)
    os.makedirs(ig_tools_init.ig_bin_directory)

def CompileAuxBinaries(log):
    log.info("\n======== Compiling binaries starts\n")
    PrepareBinDirectory()
    work = subprocess.Popen(ig_tools_init.make_ig_bins, shell = True)
    work.wait()
    log.info("\n======== Compiling binaries ends\n")

def CheckForBinExistance(log):
    if os.path.exists(ig_tools_init.PathToBins.paired_read_merger_tool):
        log.info("* Tool bin/ig_tools/paired_read_merger was successfully installed")
    else:
        log.info("ERROR: Tool bin/ig_tools/paired_read_merger was not found")
        sys.exit(1)

    if os.path.exists(ig_tools_init.PathToBins.fastq_to_fasta_tool):
        log.info("* Tool bin/ig_tools/fastq_to_fasta was successfully installed")
    else:
        log.info("ERROR: Tool bin/ig_tools/fastq_to_fasta was not found")
        sys.exit(1)

    if os.path.exists(ig_tools_init.PathToBins.merged_reads_stats_calc_tool):
        log.info("* Tool bin/ig_tools/merged_reads_stats_calc_tool was successfully installed")
    else:
        log.info("ERROR: Tool bin/ig_tools/merged_reads_stats_calc_tool was not found")
        sys.exit(1)

def RunPreparation(params, log):
    # edit config
    PrepareConfig(params, ig_tools_init.path_to_config_template, log)

    # complile aux binaries
    CompileAuxBinaries(log)

    # check for existance
    CheckForBinExistance(log)

def main():
    # preparing logging
    log = logging.getLogger('prepare_ig_tools')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    # preparing command line arguments
    try:
        options, not_options = getopt.gnu_getopt(sys.argv, BaseOptions.short_options, BaseOptions.long_options)
    except getopt.GetoptError:
        _, exc, _ = sys.exc_info()
        sys.stderr.write(str(exc) + "\n")
        usage(log)
        sys.stderr.flush()
        sys.exit(1)

    params = Params()
    for opt, arg in options:
        if opt == '--path-to-igblast':
            params.path_to_igblast = arg
        elif opt == '--help':
            usage(log)
            sys.exit(0)
        elif opt == '--dev-mode':
            params.path_to_igblast = os.path.abspath("src/tools/igblast")
    PrintParams(params, log)
    CheckForParamsCorrectness(params, log)

    # log
    log_filename = os.path.join(ig_tools_init.home_directory, "prepare_ig_tools.log")
    if os.path.exists(log_filename):
        os.remove(log_filename)
    log_handler = logging.FileHandler(log_filename, mode='a')
    log.addHandler(log_handler)
    log.info("Log will be written to " + os.path.abspath(log_filename) + "\n")

    # run all
    try:
        RunPreparation(params, log)
    except (KeyboardInterrupt):
        log.info("\nprepare_ig_tools.py was interrupted!")
    except BaseException:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            log.info("\nERROR: Exception caught. Please contact us and send " + log_filename + " file")
    
    log.info("\nLog was written to " + os.path.abspath(log_filename))

if __name__ == '__main__':
    main()
