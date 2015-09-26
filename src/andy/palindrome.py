#!/usr/bin/env python
import sys
import os
from sets import Set
from os import listdir
from os.path import isfile, isdir, join
import getopt
import logging
import shutil
import numpy

self_directory = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)))) + '/'
home_directory = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')) + '/'
repo_directory = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..')) + '/'
python_utils = os.path.join(home_directory, "ig_tools/python_utils/")

sys.path.append(python_utils)
sys.path.append(repo_directory)

import igblast_utils
import init
import pandas
    
class BaseOptions:
    long_options = "threads= max-mismatch= min-overlap= species= test skip-drawing only-merging help".split()
    short_options = "1:2:o:s:a:"

class Params:
    merged_fastq_reads = ""
    dataset_name = ""
    output_dir = ""
    igblast_align_output = ""

def PrepareOutputDirectory(output_dir_path):
    if os.path.exists(output_dir_path):
        shutil.rmtree(output_dir_path)
    os.makedirs(output_dir_path)

def FillOutputNames(params, basename):
    params.dataset_name = basename
    params.output_dir = os.path.join(self_directory, basename)

def usage(log):
    log.info("palindrome.py [options] -s merged_reads.fastq -a ig_blast.out -o <output_dir>")

    log.info("\nBasic options:")
    log.info("  -s\t\t\t<filename>\tFASTQ file with merged cleaned reads (required)")
    log.info("  -a\t\t\t<filename>\tfile containing output of IgBlast for the merged reads (-a)")
    log.info("  -o\t\t\t<directory>\toutput directory (required)")
    log.info("  --test\t\t\t\truns test dataset")
    log.info("  --help\t\t\t\tprints help")

def CheckForParams(params, log):
    if params.merged_fastq_reads == "":
        log.info("ERROR: input reads were not specified")
        usage(log)
        sys.exit(1)
    if not os.path.exists(params.merged_fastq_reads):
        log.info("ERROR: File with merged reads " + params.merged_fastq_reads + " was not found")
        usage(log)
        sys.exit(1)

def parse_command_line(options, log):
    params = Params()
    for opt, arg in options:
        if opt == '-o':
            FillOutputNames(params, arg)
        elif opt == '-s':
            params.merged_fastq_reads = os.path.abspath(arg)
        elif opt == '-a':
            params.igblast_align_output = os.path.abspath(arg)
        elif opt == '--test':
            params.merged_fastq_reads = os.path.join(repo_directory, 'test_dataset/merged_reads.fastq')
            log.info(params.merged_fastq_reads)
            params.igblast_align_output = os.path.join(repo_directory, 'test_dataset/igblast_cleaned.align')
            log.info("igblast_cleaned " + params.igblast_align_output)
            FillOutputNames(params, "palindrome_test")
        elif opt == "--help":
            usage(log)
            sys.exit(0)
    return params

def comp_let(a, b):
    str = 'ATCG'
    return str.find(a) == str.find(b) ^ 1

def GetPallindromeLen(seq, left, log):
    right = left + 1
    #log.info(str(left) + " " + str(right))
    while left >= 0 and right < len(seq) and comp_let(seq[left], seq[right]):
        left -= 1 
        right += 1
    return (right - left - 1) / 2

def SetPallindrome(stats, i, hit_row, seq, orientation, log):
    if orientation == 'right':
        clavage = hit_row.s_end != hit_row.subject_length - 1
        pal_len = GetPallindromeLen(seq, hit_row.q_end, log)
    else:
        clavage = hit_row.s_start != 0
        pal_len = GetPallindromeLen(seq, hit_row.q_start - 1, log) 
    if hit_row.type == 'D':
        type = hit_row.type + " " + orientation
    else: 
        type = hit_row.type
    stats.append({'N' : i, 'Clavage' : clavage, 'Type' : type, 'Length' : pal_len })
    #stats.append([i, clavage, type, pal_len])

def GetPallindromeStatistics(params, log):
    igblast_output = igblast_utils.ParseIgBlastOutput(params.igblast_align_output, log) 
    log.info("IgBlast output was parsed")   
    merged_fastq_lines = open(params.merged_fastq_reads, "r").readlines()
    #log.info(len(merged_fastq_lines))
    
    num_merged_reads = len(merged_fastq_lines) / 4
    #num_merged_reads = 5000
    stats = []
    log.info(str(num_merged_reads) + " merged reads were processed") 
    for i in range(0, num_merged_reads):
        name = merged_fastq_lines[i * 4].strip()[1:]
        seq = merged_fastq_lines[i * 4 + 1].strip()
        #qual = merged_fastq_lines[i * 4 + 3].strip()
        igblast_block = igblast_output[name]
        for hit_row in igblast_block.hit_table:
            #log.info(hit_row.type)
            if hit_row.type == 'V':
                SetPallindrome(stats, i, hit_row, seq, 'right', log)
            elif hit_row.type == 'D':
                SetPallindrome(stats, i, hit_row, seq, 'left', log)
                SetPallindrome(stats, i, hit_row, seq, 'right', log)
            elif hit_row.type == 'J':
                SetPallindrome(stats, i, hit_row, seq, 'left', log)
    return stats        

def main():
    # prepare log
    log = logging.getLogger('palindrome')
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
    #PrintInputParams(params, log)
    log.info("Collecting palindrome statistics started")
    stats = GetPallindromeStatistics(params, log)
    log.info("Collecting palindrome statistics finished")
    pandas.DataFrame(stats).to_csv(params.output_dir + "/statistics.csv", sep = '\t', index = False)
    log.info("Palindrome statistics is located at " + params.output_dir + "/statistics.csv")

if __name__ == '__main__':
    main()
