#!/usr/bin/env python2
import logging
import os
import sys

def CreateLogger():
    log = logging.getLogger('print_umi_graph_stats')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    return log

def ParseCommandLineParams(log):
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Prints graph decomposition statistics.")
    input_group = parser.add_mutually_exclusive_group(required=False)
    input_group.add_argument("-i", "--input", type=str, dest="input_file", help="Input file path")
    pair_reads = input_group.add_argument_group("Paired-end reads")
    pair_reads.add_argument("-1", type=str, dest="left_reads", help="Left paired-end reads in FASTQ format")
    pair_reads.add_argument("-2", type=str, dest="right_reads", help="Right paired-end reads in FASTQ format")
    parser.add_argument("-o", "--output", type=str, dest="stats_file", help="Output statistics file path", required=True)
    parser.add_argument("--tmp", type=str, dest="tmp_dir", default=".", help="Temporary files directory path")
    parser.add_argument("-c", "--clean", dest="clean", action="store_true", help="Will remove all temporary files")
    parser.add_argument("-t", "--threads", type=int, dest="threads", help="Number of threads to be used")
    parser.add_argument("--tau", type=int, dest="tau", default=1, help="Distance threshold for the UMI graph")
    parser.add_argument("--umi-cleavage-length", type=int, dest="umi_cleavage_length", default=0, help="Cleave UMIs by the specified length")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    params = parser.parse_args()

    if bool(params.left_reads) != bool(params.right_reads):
        parser.print_help()
        sys.exit(1)

    return params

class BinaryPaths:
    def __init__(self):
        run_directory = os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))) + '/'
        self.bin_dir = os.path.join(run_directory, "build/release/bin/")
        self.paired_read_merger = os.path.join(self.bin_dir, "paired_read_merger")
        self.vj_finder = os.path.join(self.bin_dir, "ig_kplus_vj_finder")
        self.umi_to_fastq = os.path.join(self.bin_dir, "umi_to_fastq")
        self.ig_trie_compressor = os.path.join(self.bin_dir, "ig_trie_compressor")
        self.ig_swgraph_construct = os.path.join(self.bin_dir, "ig_swgraph_construct")
        self.print_graph_decomposition_stats = os.path.join(self.bin_dir, "print_graph_decomposition_stats")

class BinaryRunner:
    def __init__(self, binary, input_file, output_file, params = ""):
        self.binary = binary
        self.input_file = input_file
        self.output_file = output_file
        self.params = params

    def Run(self, log, clean, threads=None):
        if not clean and os.path.exists(self.output_file):
            log.info("Output file %s already found for %s. Skipping.", self.output_file, self.binary)
            return
        if not os.path.exists(self.binary):
            log.error("Binary %s not found", self.binary)
            exit(1)
        if self.input_file and not os.path.exists(self.input_file):
            log.error("Input file %s for binary %s not found", self.input_file, self.binary)
            exit(1)
        if os.path.exists(self.output_file):
            os.remove(self.output_file)

        cmdline = "%s %s %s %s %s" % (
            self.binary,
            "-i " + self.input_file if self.input_file else "",
            "-o " + self.output_file if self.output_file else "",
            "-t " + str(threads) if threads else "",
            self.params)
        log.info("Running " + cmdline)
        exit_code = os.system(cmdline)
        if exit_code != 0:
            exit(exit_code)

class WorkflowRunner:
    def Run(self, log, params):
        # if params.clean:
        #     if os.path.exists(params.tmp_dir):
        #         log.info("Removing old tmp directory")
        #         import shutil
        #         shutil.rmtree(params.tmp_dir)
        #     log.info("Creating new tmp directory at " + params.tmp_dir)
        if not os.path.exists(params.tmp_dir):
            os.makedirs(params.tmp_dir)

        self.threads = params.threads

        if not params.input_file:
            merged_reads = self.MergeReads(log, params.left_reads, params.right_reads, params.tmp_dir, params.clean)
            single_reads = self.CleanReads(log, merged_reads, params.tmp_dir, params.clean)
        else:
            single_reads = params.input_file
        umi_fastq = self.ExtractUmi(log, single_reads, params.tmp_dir, params.umi_cleavage_length, params.clean)
        umi_compressed = self.CompressFastq(log, umi_fastq, params.tmp_dir, params.clean)
        umi_graph = self.ConstructGraph(log, umi_compressed, params.tmp_dir, params.tau, params.clean)
        self.PrintStats(log, umi_compressed, umi_graph, params.stats_file)

    def GenerateSingleReadsFileName(self, log, left_reads_file, right_reads_file, out_dir):
        left_reads = os.path.split(left_reads_file)[1]
        right_reads = os.path.split(right_reads_file)[1]
        default = os.path.join(out_dir, "single_reads" + os.path.splitext(left_reads)[1])
        if len(left_reads) != len(right_reads):
            return default
        unequal_chars = (i for i in xrange(len(left_reads)) if left_reads[i] != right_reads[i])
        index = next(unequal_chars, -1)
        if index == -1:
            return default
        if left_reads[index] != '1' or right_reads[index] != '2':
            return default
        if next(unequal_chars, -1) != -1:
            return default
        return os.path.join(out_dir, left_reads[: index] + "12" + left_reads[index + 1 :])

    def MergeReads(self, log, left_reads_file, right_reads_file, out_dir, clean):
        output_file = self.GenerateSingleReadsFileName(log, left_reads_file, right_reads_file, out_dir)
        if not clean and os.path.exists(output_file):
            log.info("Skipping reads merge.")
            return output_file
        BinaryRunner(BinaryPaths().paired_read_merger, "", "", "%s %s %s" % (left_reads_file, right_reads_file, output_file)).Run(log, clean)
        return output_file

    def CleanReads(self, log, input_file, out_dir, clean):
        result_dir = os.path.join(out_dir, "vdj_finder")
        output_file = os.path.join(result_dir, "cleaned_reads.fa")
        if not clean and os.path.exists(output_file):
            log.info("Skipping clean reads.")
            return output_file
        params = "-o " + result_dir + " --db-directory " + os.path.join(BinaryPaths().bin_dir, "germline")
        BinaryRunner(BinaryPaths().vj_finder, input_file, "", params).Run(log, clean, self.threads)
        return output_file

    def ExtractUmi(self, log, input_file, out_dir, umi_cleavage_length, clean):
        output_file = os.path.join(out_dir, os.path.splitext(os.path.split(input_file)[1])[0] + "_umi.fastq")
        BinaryRunner(BinaryPaths().umi_to_fastq, input_file, output_file, "-c %d" % umi_cleavage_length).Run(log, clean)
        return output_file

    def CompressFastq(self, log, input_file, out_dir, clean):
        output_file = os.path.join(out_dir, os.path.splitext(os.path.split(input_file)[1])[0] + "_compressed.fastq")
        BinaryRunner(BinaryPaths().ig_trie_compressor, input_file, output_file).Run(log, clean)
        return output_file

    def ConstructGraph(self, log, input_file, out_dir, tau, clean):
        output_file = os.path.join(out_dir, os.path.splitext(os.path.split(input_file)[1])[0] + ".graph")
        with open(input_file, "r") as umis:
            i = 0
            for line in umis:
                if i == 1:
                    line = line.rstrip('\n')
                    umi_length = len(line)
                    break
                i += 1
        BinaryRunner(BinaryPaths().ig_swgraph_construct, input_file, output_file, "-k %d --tau %d -A" % (umi_length // (tau + 1), tau) ).Run(log, clean, self.threads)
        return output_file

    def PrintStats(self, log, reads_file, graph_file, output_file):
        BinaryRunner(BinaryPaths().print_graph_decomposition_stats, "", output_file, "-r %s -g %s" % (reads_file, graph_file)).Run(log, True)

def main():
    log = CreateLogger()
    params = ParseCommandLineParams(log)

    WorkflowRunner().Run(log, params)

if __name__ == '__main__':
    main()
