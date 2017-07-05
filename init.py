############################################################################
# Copyright (c) 2011-2015 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import getopt
import os
import logging
import shutil
import datetime
from time import gmtime, strftime

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'
ig_bin_directory = os.path.join(home_directory, "build/release/bin/")
python_src_directory = os.path.join(home_directory, "py/utils/")
config_directory = os.path.join(home_directory, "configs/ig_tools/")
pipeline_dir = os.path.join(home_directory, "py/pipeline/")
igblast_directory = os.path.join(home_directory, "src/tools/igblast/")
ms_utils_directory = os.path.join(home_directory, "py/mass_spectra_analysis/")

path_to_config_template = os.path.join(config_directory, "config.info.template")

sys.path.append(pipeline_dir)
sys.path.append(python_src_directory)
sys.path.append(ms_utils_directory)

class PathToBins:
    paired_read_merger_tool = os.path.join(ig_bin_directory, "paired_read_merger")
    fastq_to_fasta_tool = os.path.join(ig_bin_directory, "fastq_to_fasta")
    merged_reads_stats_calc_tool = os.path.join(ig_bin_directory, "compute_merged_reads_stats")
    dense_subgraph_finder = os.path.join(ig_bin_directory, "dense_sgraph_finder")

    run_paired_read_merger_tool = ig_bin_directory + "./paired_read_merger"
    run_fastq_to_fasta_tool = ig_bin_directory + "./fastq_to_fasta"
    run_merged_reads_stats_calc_tool = ig_bin_directory + "./compute_merged_reads_stats"
    run_igblast = os.path.join(igblast_directory, "bin/igblastn")
    run_dense_sgraph_finder = os.path.join(ig_bin_directory, "./dense_sgraph_finder")

def PrintCommandLine(argv, log):
    command_line = " ".join([str(x) for x in argv] )
    log.info("Command line: "+ command_line)

def ReadConfig():
    if not os.path.exists(path_to_config_template):
        print("ERROR: config file " + path_to_config_template + " was not found")
    f = open(path_to_config_template, "r")
    config_params = dict()
    for line in f.readlines():
        splits = line.split()
        config_params[splits[0]] = splits[1]
    return config_params

def ErrorMsg(log):
    log.error("Something goes wrong. Please contact us and send .log file")
    sys.exit(1)

def AbnormalFinishMsg(log, program_name):
    log.info("Script " + program_name + " finished abnormally. Please contact us and send .log file")
