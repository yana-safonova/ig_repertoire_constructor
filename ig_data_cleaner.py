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

import init
import drawing_utils
import igblast_utils
import files_utils

config = {}

class BaseOptions:
    long_options = "threads= max-mismatch= min-overlap= species= test skip-drawing only-merging help".split()
    short_options = "1:2:o:s:a:"

class Params:
    left_reads = ""
    right_reads = ""
    merged_fastq_reads = ""
    merged_fasta_reads = ""
    max_mismatch = 0.1
    min_overlap = 50
    min_length = 300
    species = 'human'
    dataset_name = ""
    output_dir = ""
    basename = ""
    draw_stats = True
    only_merging = False

    start_from_merging = True
    start_from_cleaning = False

    igblast_align_output = ""
    bad_reads = Set()

    cleaned_reads = ""
    filtered_reads = ""
    igblast_filtered_align_output = ""

    aver_paired_qual_stats = ""
    aver_merged_qual_stats = ""
    merged_rl_stats = ""

    rl_hist = ""
    aver_quality_plot = ""

    log = ""

def PrintInputParams(params, log):
    log.info("\nInput parameters:")
    log.info("Output directory:\t\t\t" + params.output_dir)
    if params.start_from_merging:
        log.info("Left reads:\t\t\t\t" + params.left_reads)
        log.info("Right reads:\t\t\t\t" + params.right_reads)
    if params.start_from_cleaning:
        log.info("Merged reads:\t\t\t\t" + params.merged_fastq_reads)
    if params.igblast_align_output != "":
        log.info("IgBlast output:\t\t\t" + params.igblast_align_output)
    log.info("Min allowed overlap size:\t\t" + str(params.min_overlap))
    log.info("Min allowed read length:\t\t" + str(params.min_length))
    log.info("Max allowed mismatch rate:\t\t" + str(params.max_mismatch))
    log.info("Species:\t\t\t\t" + params.species)
    if params.draw_stats:
        log.info("Drawing plots is enable")
    else:
        log.info("Drawing plots is disable")

def usage(log):
    log.info("ig_data_cleaner.py [options] -1 left_reads.fastq -2 right_reads.fastq -o <output_dir>")
    log.info("OR:")
    log.info("ig_data_cleaner.py [options] -s merged_reads.fastq -o <output_dir>")

    log.info("\nBasic options:")
    log.info("  -1\t\t\t<filename>\tFASTQ file with left reads (required)")
    log.info("  -2\t\t\t<filename>\tFASTQ file with right reads (required)")
    log.info("  -s\t\t\t<filename>\tFASTQ file with merged uncleaned reads (required)")
    log.info("  -o\t\t\t<directory>\toutput directory (required)")
    log.info("  --test\t\t\t\truns test dataset")

    log.info("\nAdvanced options:")
    log.info("  -a\t\t\t<filename>\tfile containing output of IgBlast for the merged reads (-a)")
    log.info("  --min-overlap\t\t<int>\t\tminimal allowed size of overlap in paired reads merging [default: '60']")
    log.info("  --min-length\t\t<int>\t\tminimal allowed length of merged read [default: '300']")
    log.info("  --max-mismatch\t<float>\t\tmaximal allowed mismatch of overlap in paired reads merging [default: '0.1']")
    log.info("  --species\t\t\t\tspecies, possible values: human, mouse [default: 'human']")
    log.info("  --only-merging\t\t\truns only paired read merger and visualization of stats")
    log.info("  --skip-drawing\t\t\tskips visualization of statistics for merged reads")
    log.info("  --help\t\t\t\tprints help")

def CheckForParams(params, log):
    if params.left_reads == "" and params.right_reads == "" and params.merged_fastq_reads == "":
        log.info("ERROR: input reads were not specified")
        usage(log)
        sys.exit(1)
    else:
        if params.left_reads == "" and params.right_reads != "":
            log.info("ERROR: left reads (-1) were not specified")
            usage(log)
            sys.exit(1)
        if params.left_reads != "" and params.right_reads == "":
            log.info("ERROR: right reads (-2) were not specified")
            usage(log)
            sys.exit(1)
        if params.left_reads != "" and params.right_reads != "" and params.merged_fastq_reads != "":
            log.info("WARN: both paired-end and merged reads were specified. ig_data_preparation will use " + params.left_reads + " and " + params.right_reads)
    if params.start_from_merging:
        if not os.path.exists(params.left_reads):
            log.info("ERROR: File with left reads " + params.left_reads + " was not found")
            usage(log)
            sys.exit(1)
        if not os.path.exists(params.right_reads):
            log.info("ERROR: File with right reads " + params.right_reads + " was not found")
            usage(log)
            sys.exit(1)
    if params.start_from_cleaning:
        if not os.path.exists(params.merged_fastq_reads):
            log.info("ERROR: File with merged reads " + params.merged_fastq_reads + " was not found")
            usage(log)
            sys.exit(1)
    if params.max_mismatch < 0 or params.max_mismatch > 1:
        log.info("ERROR: Maximal allowed mismatch rate (--max-mismatch) should be from [0, 1]")
        usage(log)
        sys.exit(1)
    if params.species not in ['human', 'mouse']:
        log.info('ERROR: Species must be human or mouse')
        usage(log)
        sys.exit(1)

def FillOutputNames(params, basename):
    params.dataset_name = basename
    params.output_dir = os.path.join(ig_tools_init.home_directory, basename)
    params.basename = os.path.join(params.output_dir, params.dataset_name)
    
def parse_command_line(options, log):
    params = Params()
    for opt, arg in options:
        if opt == '-o':
            FillOutputNames(params, arg)
        elif opt == '-1':
            params.left_reads = os.path.abspath(arg)
        elif opt == '-2':
            params.right_reads = os.path.abspath(arg)
        elif opt == '-s':
            params.merged_fastq_reads = os.path.abspath(arg)
            params.start_from_merging = False
            params.start_from_cleaning = True
        elif opt == '-a':
            params.igblast_align_output = os.path.abspath(arg)
        elif opt == "--max-mismatch":
            params.max_mismatch = int(arg)
        elif opt == "--min-overlap":
            params.min_overlap = int(arg)
        elif opt == "--min-length":
            params.min_length = int(arg)
        elif opt == '--species':
            params.species = arg
        elif opt == '--skip-drawing':
            params.draw_stats = False
        elif opt == '--test':
            params.left_reads = os.path.join(ig_tools_init.home_directory, 'ig_test_dataset/left.fastq')
            params.right_reads = os.path.join(ig_tools_init.home_directory, 'ig_test_dataset/right.fastq')
            FillOutputNames(params, "ig_data_cleaner_test")
        elif opt == '--only-merging':
            params.only_merging = True
        elif opt == "--help":
            usage(log)
            sys.exit(0)
    return params

def PrepareOutputDirectory(output_dir_path):
    if os.path.exists(output_dir_path):
        shutil.rmtree(output_dir_path)
    os.makedirs(output_dir_path)

# --------------------------- PairedReadMerger --------------------------

def MergePairedReads(params, log):
    if params.start_from_cleaning:
        return 
    log.info("==== Merging paired-end reads")
    command_line = ig_tools_init.PathToBins.run_paired_read_merger_tool + " " + params.left_reads + " " + params.right_reads + " " + params.output_dir + "/merged_reads" + " --max-mismatch=" + str(params.max_mismatch) + "  --min-overlap=" + str(params.min_overlap)
    error_code = os.system(command_line + " 2>&1 | tee -a " + params.log)

    if error_code != 0:
        AbnormalFinishMsg(log, "paired_read_merged")
        sys.exit(1)

    params.merged_fastq_reads = os.path.join(params.output_dir, "merged_reads.fastq")
    if os.path.exists(params.merged_fastq_reads):
        log.info("* Merged reads were written to " + params.merged_fastq_reads)
    else:
        log.info("ERROR: FASTQ with merged reads was not found")
        ig_tools_init.ErrorMsg(log)

# --------------------------- ContaminationCleaning --------------------------

def GetFastaName(fastq_name):
    splits = fastq_name.split('.')
    fasta_name = ""
    for i in range(len(splits) - 1):
         fasta_name += splits[i] + "."
    return fasta_name + "fasta"

def FastqToFasta(params, log):
    log.info("\n==== Conversion from FASTQ to FASTA")
    params.merged_fasta_reads = GetFastaName(params.merged_fastq_reads)
    error_code = os.system(ig_tools_init.PathToBins.run_fastq_to_fasta_tool + " " + params.merged_fastq_reads + " " + params.merged_fasta_reads + " 2>&1 | tee -a " + params.log)

    if error_code != 0:
        AbnormalFinishMsg(log, "fastq_to_fasta")
        sys.exit(1)

    if os.path.exists(params.merged_fasta_reads):
        log.info("* FASTA file with merged reads was written to " + params.merged_fasta_reads)
    else:
        log.info("FASTA file with merged reads was not found")
        ig_tools_init.ErrorMsg(log)

def RunIgblast(params, log):
    if params.igblast_align_output != "":
        return

    log.info("\n==== Running IgBLAST")
    params.igblast_align_output = os.path.join(params.output_dir, "igblast.align")
    igblast_command_line = ig_tools_init.RunIgblast() + " -germline_db_V "+ ig_tools_init.IgblastDirectory() +"database/" + params.species + "_gl_V -germline_db_J "+ ig_tools_init.IgblastDirectory() +"database/" + params.species + "_gl_J -germline_db_D "+ ig_tools_init.IgblastDirectory() +\
             "database/" + params.species + "_gl_D -query "+ params.merged_fasta_reads + " -show_translation -auxiliary_data auxilary_file -num_alignments 10  -outfmt \"7 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore slen\" -domain_system imgt > "\
            + params.igblast_align_output
    os.environ['IGDATA'] = ig_tools_init.IgblastDirectory()
    error_code = os.system(igblast_command_line + " 2>&1 | tee -a " + params.log)

    if error_code != 0:
        AbnormalFinishMsg(log, "IgBlast")
        sys.exit(1)

    if os.path.exists(params.igblast_align_output):
        log.info("* Output of IgBLAST alignment was written to " + params.igblast_align_output)
    else:
        log.info("ERROR: Output of IgBLAST alignment was not found")
        ig_tools_init.ErrorMsg(log)
    
def FilterBadReads(params, log):
    log.info("\n==== Searching for contaminated reads")
    bf = open(params.igblast_align_output, "r")
    curr_read = ""
    line = "smth"
    while line:
        while curr_read == "" :
            line = bf.readline()
            if not line:
                break
            arr = line.split()
            if len(arr) >=3 and arr[1] == "Query:":
                curr_read = arr[2]
                break
        if not line:
            break
        while line and line.find("hits found") == -1:
            line = bf.readline()
        hitnum = int(line.split(' ')[1])
        aligned = False
        for i in range(0, hitnum ):
            line = bf.readline()
            e_value = float(line.split()[12])
            if e_value < 0.001:
                aligned = True
        if not aligned:
            params.bad_reads.add(curr_read)
        curr_read = ""
    log.info(str(len(params.bad_reads)) + " read(s) will be filtered")

def FilterShortReads(params, log):
    log.info("\n=== Searching for short reads")
    merged_fastq_file = open(params.merged_fastq_reads, "r")
    curr_read = ""
    count = 0
    bad_reads_prev_count = len(params.bad_reads)
    for line in merged_fastq_file:
        if count % 4 == 0:
            curr_read = line.strip()[1:]
        if count % 4 == 1 and len(line.strip()) < params.min_length:
            params.bad_reads.add(curr_read)
        count += 1
    log.info(str(len(params.bad_reads) - bad_reads_prev_count) + " read(s) will be filtered")

def WriteCleanedFilteredReads(params, log):
    log.info("\n==== Writing cleaned and filtered reads")
    params.cleaned_reads = os.path.join(params.output_dir, "cleaned_reads.fastq")
    params.filtered_reads = os.path.join(params.output_dir, "filtered_reads.fastq")
    merged_fastq_lines = open(params.merged_fastq_reads, "r").readlines()
    igblast_output = igblast_utils.ParseIgBlastOutput(params.igblast_align_output, log) 
    log.info("IgBlast output was parsed")   

    num_filtered = 0
    num_cleaned = 0

    num_merged_reads = len(merged_fastq_lines) / 4
    for i in range(0, num_merged_reads):
        name = merged_fastq_lines[i * 4].strip()[1:]
        seq = merged_fastq_lines[i * 4 + 1].strip()
        qual = merged_fastq_lines[i * 4 + 3].strip()

        reverse = not igblast_output.GetBlockByName(name).vdj_rearrangement.direct_strand
#        print "Read: " + name + ", Reverse: " + str(reverse)

        if name in params.bad_reads:
            files_utils.WriteReadInFastqFile(name, seq, qual, reverse, params.filtered_reads)
            num_filtered += 1
        else:
            files_utils.WriteReadInFastqFile(name, seq, qual, reverse, params.cleaned_reads)
            num_cleaned += 1

    if os.path.exists(params.cleaned_reads):
        log.info("* " + str(num_cleaned) + " cleaned reads were written to " + params.cleaned_reads)
    else:
        log.info("ERROR: cleaned reads were not found")
        ig_tools_init.ErrorMsg(log)

    if os.path.exists(params.filtered_reads):
        log.info("* " + str(num_filtered) + " filtered reads were written to " + params.filtered_reads)    
    else:
        log.info("ERROR: filtered reads were not found")
        ig_tools_init.ErrorMsg(log)


def WriteFilteredAlignOutput(params, log):
    log.info("\n==== Correction of IgBLAST alignment output for cleaned reads")
    count = 0
    params.igblast_cleaned_align_output = os.path.join(params.output_dir, "igblast_cleaned.align")
    all_align_output = open(params.igblast_align_output, "r")
    cleaned_align_output = open(params.igblast_cleaned_align_output, "w")
    curr_read = ""
    for line in all_align_output:
        arr = line.split()
        if len(arr) >= 3 and arr[1] == "Query:":
            curr_read = arr[2]
        if not curr_read in params.bad_reads:
            cleaned_align_output.write(line)
    log.info("* IgBLAST align output for cleaned reads was written to " + params.igblast_cleaned_align_output)    

def RunContaminationCleaning(params, log):
    if params.only_merging:
        params.cleaned_reads = params.merged_fastq_reads
        params.filtered_reads = params.merged_fastq_reads
        return 

    FastqToFasta(params, log)
    RunIgblast(params, log)
    FilterBadReads(params, log)
    FilterShortReads(params, log)
    WriteCleanedFilteredReads(params, log)
    WriteFilteredAlignOutput(params, log)

# --------------------------- compute stats -----------------

def aver(data_list):
    aver = 0
    for item in data_list:
        aver += item / len(data_list)
    return aver

def ProcessAverageQualityStats(params, log):
    paired_read_data = drawing_utils.ReadFloatGraphicalData(params.aver_paired_qual_stats)
    merged_read_data = drawing_utils.ReadFloatGraphicalData(params.aver_merged_qual_stats)
    log.info("Average quality of original paired reads\t%.2f" %aver(paired_read_data.all_keys))
    log.info("Average quality of merged reads\t\t\t%.2f" %aver(merged_read_data.all_keys))
    log.info("Average quality improvement\t\t\t%.2f" %(aver(merged_read_data.all_keys) - aver(paired_read_data.all_keys)) + "\n")
    x = list()
    x.append(range(0, len(paired_read_data.all_keys)))
    x.append(range(0, len(merged_read_data.all_keys)))
    y = list()
    y.append(paired_read_data.all_keys)
    y.append(merged_read_data.all_keys)

    params.aver_qual_plot = os.path.join(params.output_dir, "aver_quality.png")
    setting = drawing_utils.GetGraphicalSettings(xlabel = "Nucleotide position", ylabel = "Average quality", title = "", output_filename = params.aver_qual_plot, bins = 0, label = ["Original paired reads", "Merged reads"], histtype = None, xlog_scale = False, ylog_scale = False, draw_legend = True, colors = ["r", 'b'], legend_loc = 'lower right') 
    drawing_utils.DrawMultiplePlot(x, y, setting)
    if os.path.exists(params.aver_qual_plot):
        log.info("* Graphics of average nucleotide quality for original paired and merged reads were written to " + params.aver_qual_plot)    
    else:
        log.info("ERROR: Graphics of average nucleotide quality for original paired and merged reads were not found")
        ig_tools_init.ErrorMsg(log)

    
def DrawReadLengthHistogram(params, log):
    hist = drawing_utils.ReadIntGraphicalData(params.merged_rl_stats)
    num_bars = int((hist.max_cluster - hist.min_cluster) * .75) + 1

    params.rl_hist = os.path.join(params.output_dir, "merged_read_length.png")
    hist_setting = drawing_utils.GetGraphicalSettings(xlabel = 'Read length', ylabel = 'Read number', title = "", output_filename = params.rl_hist, bins = num_bars)  
    drawing_utils.DrawHistogram(hist.all_keys, hist_setting)

    if os.path.exists(params.rl_hist):
        log.info("* Histogram of merged read length distribution was written to " + params.rl_hist + "\n")
    else:
        log.info("ERROR: Histogram of merged read length distribution was not found")
        ig_tools_init.ErrorMsg(log)

def RunComputeStats(params, log):
    if not params.start_from_merging:
        return 
    log.info("\n==== Calculatiom of statistics for merged reads")
    command_line = ig_tools_init.PathToBins.run_merged_reads_stats_calc_tool + " " + params.left_reads + " " + params.right_reads + " " + params.cleaned_reads + " " + params.output_dir
    error_code = os.system(command_line + " 2>&1 | tee -a " + params.log)

    if error_code != 0:
        AbnormalFinishMsg(log, "merged_reads_stats_calculator")
        sys.exit(1)

    params.aver_merged_qual_stats = os.path.join(params.output_dir, "merged_nucl_qual.stats")
    if os.path.exists(params.aver_merged_qual_stats):
        log.info("* Statistics about average merged reads quality were written to " + params.aver_merged_qual_stats)
    else:
        log.info("ERROR: Statistics about average merged reads quality were not found")
        ig_tools_init.ErrorMsg(log)

    params.aver_paired_qual_stats = os.path.join(params.output_dir, "paired_nucl_qual.stats")
    if os.path.exists(params.aver_paired_qual_stats):
        log.info("* Statistics about average paired reads quality were written to " + params.aver_paired_qual_stats)
    else:
        log.info("ERROR: Statistics about average paired reads quality were not found")
        ig_tools_init.ErrorMsg(log)

    params.merged_rl_stats = os.path.join(params.output_dir, "merged_rl.stats")
    if os.path.exists(params.merged_rl_stats):
        log.info("* Statistics about merged read lengths were written to " + params.merged_rl_stats)
    else:
        log.info("ERROR: Statistics about merged read lengths were not found")
        ig_tools_init.ErrorMsg(log)

    DrawReadLengthHistogram(params, log)
    ProcessAverageQualityStats(params, log)
   
# --------------------------- final output ------------------
def CleanOutputDir(params):
    if os.path.exists(params.merged_fasta_reads):
        os.remove(params.merged_fasta_reads)

def FinalOutput(params, log):
    log.info("Main output files:")
    log.info("* Cleaned merged reads were written to " + params.cleaned_reads)
    if params.igblast_filtered_align_output != "":
        log.info("* IgBlast alignment output for cleaned reads was written to " + params.igblast_cleaned_align_output)
    log.info("\nThank you for using IgDataCleaner!")

# --------------------------- main --------------------------

def RunIgDataCleaner(params, log):
    log.info("\n======== IgDataCleaner starts\n")
    # run paired read merger
    MergePairedReads(params, log)
    # run igblast
    RunContaminationCleaning(params, log)
    # run stats calculator
    RunComputeStats(params, log)
    log.info("\n======== IgDataCleaner ends\n")
    CleanOutputDir(params)
    FinalOutput(params, log)

def main():
    # prepare log
    log = logging.getLogger('ig_data_cleaner')
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

    # parse command line
    params = parse_command_line(options, log)

    # create output directory    
    if params.output_dir == "":
        log.info("ERROR: output directory was not specified")
        usage(log)
        sys.exit(1)
    PrepareOutputDirectory(params.output_dir)    

    log_filename = os.path.join(params.output_dir, "ig_data_cleaner.log")
    if os.path.exists(log_filename):
        os.remove(log_filename)
    log_handler = logging.FileHandler(log_filename, mode='a')
    log.addHandler(log_handler)
    params.log = log_filename
    log.info("Log will be written to " + log_filename + "\n")

    CheckForParams(params, log)

    # print input params
    ig_tools_init.PrintCommandLine(sys.argv, log)
    PrintInputParams(params, log)

    try:
        RunIgDataCleaner(params, log) 
    except (KeyboardInterrupt):
        log.info("\nIgDataCleaner was interrupted!")
    except BaseException:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            log.info("\nERROR: Exception caught. Please contact us and send .log file")

    log.info("\nLog was written to " + log_filename)

if __name__ == '__main__':
    main()
