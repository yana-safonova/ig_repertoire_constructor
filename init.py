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
import datetime
from time import gmtime, strftime

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'
ig_bin_directory = os.path.join(home_directory, "bin/ig_tools/")
python_src_directory = os.path.join(home_directory, "src/ig_tools/python_utils/")
config_directory = os.path.join(home_directory, "configs/ig_tools/")
spades_py_scripts_directory = os.path.join(home_directory, "src/spades_pipeline/")
spades_py_scripts_directory = os.path.join(home_directory, "src/spades_pipeline/")

path_to_config_template = os.path.join(config_directory, "config.info.template")

make_ig_bins = "make igtools"

sys.path.append(spades_py_scripts_directory)
sys.path.append(python_src_directory)

class PathToBins:
    paired_read_merger_tool = os.path.join(ig_bin_directory, "paired_read_merger")
    fastq_to_fasta_tool = os.path.join(ig_bin_directory, "fastq_to_fasta")
    merged_reads_stats_calc_tool = os.path.join(ig_bin_directory, "compute_merged_reads_stats")

    run_paired_read_merger_tool = ig_bin_directory + "./paired_read_merger"
    run_fastq_to_fasta_tool = ig_bin_directory + "./fastq_to_fasta"
    run_merged_reads_stats_calc_tool = ig_bin_directory + "./compute_merged_reads_stats"

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

def IgblastDirectory():
    config_params = ReadConfig()
    return config_params['path_to_igblast'] + "/"

def RunIgblast():
    config_params = ReadConfig()
    return config_params['path_to_igblast'] + "/bin/igblastn"

def ErrorMsg(log):
    log.error("Something goes wrong. Please contact us and send .log file")
    sys.exit(1)

def AbnormalFinishMsg(log, program_name):
    log.info("Script " + program_name + " finished abnormally. Please contact us and send .log file")
