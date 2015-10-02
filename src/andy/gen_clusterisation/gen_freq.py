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
from Bio import SeqIO

self_directory = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)))) + '/'
user_directory = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')) + '/'
home_directory = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..')) + '/'
repo_directory = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', '..')) + '/'
python_utils = os.path.join(home_directory, "ig_tools/python_utils/")

sys.path.append(python_utils)
sys.path.append(repo_directory)

import igblast_utils
import init
import pandas
import copy
    
class BaseOptions:
    long_options = "threads= max-mismatch= min-overlap= species= test skip-drawing only-merging help".split()
    short_options = "1:2:o:s:a:f:"

class Params:
    formatting = ""
    merged_fastq_reads = ""
    dataset_name = ""
    output_dir = ""
    output_file = ""
    igblast_align_output = ""

def PrepareOutputDirectory(output_dir_path):
    if os.path.exists(output_dir_path):
    #shutil.rmtree(output_dir_path)
        return
    os.makedirs(output_dir_path)

#def FillOutputNames(params, basename):
#    params.dataset_name = basename
#    params.output_dir = os.path.join(self_directory, basename)

def usage(log):
    log.info("palindrome.py [options] -s merged_reads.fastq -a ig_blast.out -o <output_dir>")

    log.info("\nBasic options:")
    log.info("  -f\t\t\t<format_name> FASTA or FASTQ (required)")
    log.info("  -s\t\t\t<filename>\tFASTQ file with merged cleaned reads (required)")
    log.info("  -a\t\t\t<filename>\tfile containing output of IgBlast for the merged reads (-a)")
    log.info("  -o\t\t\t<directory>\toutput directory (required)")
    log.info("  --test\t\t\t\truns test dataset")
    log.info("  --help\t\t\t\tprints help")

def CheckForParams(params, log):
    if params.formatting.lower() != "fasta" and params.formatting.lower() != "fastq":
        log.info("ERROR: input format was not specified : " + str(params.formatting.lower()))
        usage(log)
        sys.exit(1)
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
        if opt == "-f":
            #print "f " + arg
            params.formatting = arg.lower()
        elif opt == '-o':
            params.output_dir, params.output_file = os.path.split(arg)
        elif opt == '-s':
            filename, file_extension = os.path.splitext(arg)
            params.formatting = file_extension.strip()[1:]
            if params.formatting == "fa":
                params.formatting = "fasta"
            params.merged_fastq_reads = os.path.abspath(arg)
        elif opt == '-a':
            params.igblast_align_output = os.path.abspath(arg)
        elif opt == '--test':
            params.merged_fastq_reads = os.path.join(repo_directory, 'test_dataset/merged_reads.fastq')
            log.info(params.merged_fastq_reads)
            params.igblast_align_output = os.path.join(repo_directory, 'test_dataset/igblast_cleaned.align')
            log.info("igblast_cleaned " + params.igblast_align_output)
            params.output_dir = os.path.join(user_directory, 'gen_freq')
            params.output_file = 'gen_freq'
            params.formatting = 'fastq'
        elif opt == "--help":
            usage(log)
            sys.exit(0)
    return params

def comp_let(a, b):
    str = 'ATCG'
    return str.find(a) == str.find(b) ^ 1

def GetPalindromeLen(seq, left):
    right = left + 1
    #log.info(str(left) + " " + str(right))
    while left >= 0 and right < len(seq) and comp_let(seq[left], seq[right]):
        left -= 1 
        right += 1
    return (right - left - 1) / 2

def GenFrequences(params, log):
    igblast_output = igblast_utils.ParseIgBlastOutput(params.igblast_align_output, log) 
    log.info("IgBlast output was parsed")   
    #merged_fastq_lines = open(params.merged_fastq_reads, "r").readlines()
    merged_fastq_records  = list(SeqIO.parse(params.merged_fastq_reads, params.formatting))
    #log.info(len(merged_fastq_lines))
   
    #num_merged_reads = 5000
    stats = []
    for k, record in enumerate(merged_fastq_records):
        name = record.name
        #abundance = int(name[name.find('size') + 7:])
        #print name
        #abundance = int(name[name.rfind(":") + 1 : name.rfind("_")])
        abundance = int(name[name.rfind("_") + 1:])
        #print abundance
        #seq = merged_fastq_lines[i * len_format + 1].strip() # for fastq
        seq = record.seq
        prev_hit_row_type = ''
        igblast_block = igblast_output[name]
        for hit_row in igblast_block.hit_table:
            if prev_hit_row_type == hit_row.type:
                continue
            subject_id = hit_row.subject_id
            if len(subject_id) >= 3 and subject_id[-3] == '*':
                subject_id = subject_id[:-3]
            stat = {'N' : k, 'Type' : hit_row.type, 'subject_id' : subject_id, 'abundance' : abundance, 'gen_length' : hit_row.subject_length, 'perc_identity' : hit_row.perc_identity}

            if hit_row.type == 'V':
                clavage = hit_row.s_end != hit_row.subject_length - 1
                stat['Clavage'] = clavage
                stat['palindrome_len'] = GetPalindromeLen(seq, hit_row.q_end)
                stats.append(stat)
            elif hit_row.type == 'D':
                clavage = hit_row.s_start != 0
                stat['Clavage'] = clavage
                stat['Type'] = "D left" 
                stat['palindrome_len'] = GetPalindromeLen(seq, hit_row.q_start - 1)
                stats.append(stat.copy())
                clavage = hit_row.s_end != hit_row.subject_length - 1
                stat['Clavage'] = clavage
                stat['Type'] = "D right" 
                stat['palindrome_len'] = GetPalindromeLen(seq, hit_row.q_end)
                stats.append(stat)
            else:
                clavage = hit_row.s_start != 0
                stat['Clavage'] = clavage
                stat['palindrome_len'] = GetPalindromeLen(seq, hit_row.q_start - 1)
                stats.append(stat)
            prev_hit_row_type = hit_row.type
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
    log.info("Collecting Gens Frequencies started")
    stats = GenFrequences(params, log)
    log.info("Collecting Gens Frequencies finished")
    #print >> stats['V']
    #for x in stats['V']:
    #    print (x, ' ', stats[x])
    pandas.DataFrame(stats).to_csv(os.path.join(params.output_dir, params.output_file) + ".csv", sep = '\t', index = False)
    log.info("Palindrome statistics is located at " + params.output_dir + "/gen_freq.csv")


if __name__ == '__main__':
    main()
