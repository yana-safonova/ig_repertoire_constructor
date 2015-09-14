#!/usr/bin/env python2

import sys
import os
from sets import Set
from os import listdir
from os.path import isfile, isdir, join
import getopt
import logging
import shutil
import numpy
import math
from Bio import SeqIO

repo_directory = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', '..')) + '/'
user_directory = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')) + '/'
src_directory = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..')) + '/'
python_utils = os.path.join(src_directory, "ig_tools/python_utils/")

sys.path.append(repo_directory)

import init

class BaseOptions:
    long_options = "threads= max-mismatch= min-overlap= species= test skip-drawing only-merging help".split()
    short_options = "1:2:f:i:o:n:"

class Params:
    formatting = ""
    merged_reads = ""
    output_dir = ""
    n_output_files = 8 

def PrepareOutputDirectory(output_dir_path):
    if os.path.exists(output_dir_path):
        #shutil.rmtree(output_dir_path)
        return
    os.makedirs(output_dir_path)
    if os.path.exists(output_dir_path):
        shutil.rmtree(output_dir_path)
    os.makedirs(output_dir_path)

def usage(log):
    log.info("cut_merged_reads.py -f <format_name> -i <intput_dir> -o <output_dir> -n <number_of_output_files>")

    log.info("\nBasic options:")
    log.info("  -f\t\t\t<format_name> FASTA or FASTQ (required)")
    log.info("  -i \t\t\t<filename>\tFASTA/FASTQ file with merged cleaned reads")
    log.info("  -o\t\t\t<filename>\toutput filename (required)")
    log.info("  -n\t\t\t<integer>\t number of output files (default = 8)")
    log.info("  --test\t\t\t\truns test dataset")
    log.info("  --help\t\t\t\tprints help")

def CheckForParams(params, log):
    if params.formatting.lower() != "fasta" and params.formatting.lower() != "fastq":
        log.info("ERROR: input format was not specified")
        usage(log)
        sys.exit(1)
    if not os.path.exists(params.merged_reads):
        log.info("ERROR: File with merged reads " + params.merged_reads + " was not found")
        usage(log)
        sys.exit(1)

def parse_command_line(options, log):
    params = Params()
    for opt, arg in options:
        if opt == "-f":
            print "f " + arg
            params.formatting = arg.lower()
        elif opt == "-i":
            print "i " + arg
            params.merged_reads = os.path.abspath(arg)
            filename, file_extension = os.path.splitext(arg)
            params.formatting = file_extension.strip()[1:]
            #log.info("FORMAT: " + params.formatting)
        elif opt == "-o":
            print "o " + arg
            params.output_dir = os.path.abspath(arg)
        elif opt == "-n":
            print "n " + arg
            params.n_output_files = int(arg)
        elif opt == "--test":
            params.formatting = "fastq"
            params.merged_reads = os.path.join(repo_directory, 'test_dataset/merged_reads.fastq')
            log.info(params.merged_reads)
            params.output_dir = os.path.join(repo_directory, 'test_dataset/merged_reads_cut')
            #log.info("Output_directory " + params.output_dir)
            #params.n_output_files = 8
        elif opt == "--help":
            usage(log)
            sys.exit(0)
    return params

def CutReads(params, log):
    records = list(SeqIO.parse(os.path.abspath(params.merged_reads), params.formatting))
    n_records_per_file = int(math.ceil(len(records) * 1. / params.n_output_files))
    output_files = ["%s/merged_reads_cut_%.3d.%s"%(params.output_dir, j, params.formatting) 
            for j in xrange(0, params.n_output_files)]
    for i in xrange(0, len(records), n_records_per_file):
        SeqIO.write(records[i : i + n_records_per_file], 
                "%s/merged_reads_cut_%.3d.fasta"%(params.output_dir, i / n_records_per_file),
                "fasta")
    #log.info(" ".join(output_files))

def main():
    #print "repo_directory " + repo_directory + "\n"
    #print "src_directory " + src_directory + "\n"
    # prepare log
    log = logging.getLogger('cut_merged_reads')
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
    if not options:
        usage(log)
        sys.stderr.flush()
        sys.exit(1)
    
    params = parse_command_line(options, log)

    # create output directory    
    if params.output_dir == "":
        log.info("ERROR: output directory was not specified")
        usage(log)
        sys.exit(1)
    PrepareOutputDirectory(params.output_dir)    

    log_filename = os.path.join(params.output_dir, "palindrome.log")
    if os.path.exists(log_filename):
        os.remove(log_filename)
    log_handler = logging.FileHandler(log_filename, mode='a')
    log.addHandler(log_handler)
    params.log = log_filename
    log.info("Log will be written to " + log_filename + "\n")

    CheckForParams(params, log)

    # print input params
    init.PrintCommandLine(sys.argv, log)
    log.info("Cutting merged_reads...")
    CutReads(params, log)
    log.info("Cutting done. Output files are located at " + params.output_dir)

if __name__ == '__main__':
    main()
