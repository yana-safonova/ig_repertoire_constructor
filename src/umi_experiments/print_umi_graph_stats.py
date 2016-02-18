#!/usr/bin/env python2

#  ./ig_repertoire_constructor/build/release/bin/umi_to_fastq /Johnny/data/input/Ig/ibh_datasets/age_datasets/raw_reads/merged_reads/age_ig_s1_R12.fastq data/age/age_ig_s1_R12_umi.fastq
# 1611497 barcodes were extracted from /Johnny/data/input/Ig/ibh_datasets/age_datasets/raw_reads/merged_reads/age_ig_s1_R12.fastq to data/age/age_ig_s1_R12_umi.fastq
#
# sbankevich@ace:~$ ./ig_repertoire_constructor/build/release/bin/ig_trie_compressor -i data/age/age_ig_s1_R12_umi.fastq -o data/age/age_ig_s1_R12_umi_compressed.fastq
#
# sbankevich@ace:~$ ./ig_repertoire_constructor/build/release/bin/ig_swgraph_construct -k 6 --tau 1 -A -i data/age/age_ig_s1_R12_umi_compressed.fastq -o data/age/age_ig_s1_R12_umi_compressed.graph
#
# sbankevich@ace:~$ ./ig_repertoire_constructor/build/release/bin/print_graph_decomposition_stats data/age/age_ig_s1_R12_umi_compressed.graph 20
import logging
import os
import sys

run_directory = os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))) + '/'

def CreateLogger():
    log = logging.getLogger('ig_repertoire_constructor')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    return log

def ParseCommandLineParams(log):
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Prints graph decomposition statistics.")
    parser.add_argument("-i", "--input", type=str, dest="input_file", help="Input file path", required=True)
    parser.add_argument("-o", "--output", type=str, dest="stats_file", help="Output statistics file path", required=True)
    parser.add_argument("-t", "--tmp", type=str, dest="tmp_dir", default=".", help="Temporary files directory path")

    params = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    return params

class BinaryPaths:
    def __init__(self):
        bin_dir = "build/release/bin/"
        self.umi_to_fastq = os.path.join(run_directory, bin_dir, "umi_to_fastq")
        self.ig_trie_compressor = os.path.join(run_directory, bin_dir, "ig_trie_compressor")
        self.ig_swgraph_construct = os.path.join(run_directory, bin_dir, "ig_swgraph_construct")
        self.print_graph_decomposition_stats = os.path.join(run_directory, bin_dir, "print_graph_decomposition_stats")

class BinaryRunner:
    def __init__(self, binary, input_file, output_file, params = ""):
        self.binary = binary
        self.input_file = input_file
        self.output_file = output_file
        self.params = params

    def Run(self, log):
        if not os.path.exists(self.binary):
            log.error("Binary %s not found", self.binary)
            exit(1)
        if self.input_file and not os.path.exists(self.input_file):
            log.error("Input file %s for binary %s not found", self.binary, self.input_file)
            exit(1)
        if os.path.exists(self.output_file):
            os.remove(self.output_file)

        cmdline = "%s %s%s -o %s %s" % (self.binary, "-i " if self.input_file else "", self.input_file, self.output_file, self.params)
        os.system(cmdline)

class WorkflowRunner:
    def Run(self, log, params):
        if (os.path.exists(params.tmp_dir)):
            log.info("Removing old tmp directory")
            import shutil
            shutil.rmtree(params.tmp_dir)
        log.info("Creating new tmp directory at " + params.tmp_dir)
        os.makedirs(params.tmp_dir)
        umi_fastq = self.ExtractUmi(log, params.input_file, params.tmp_dir)
        umi_compressed = self.CompressFastq(log, umi_fastq, params.tmp_dir)
        umi_graph = self.ConstructGraph(log, umi_compressed, params.tmp_dir)
        self.PrintStats(log, umi_compressed, umi_graph, params.stats_file)

    def ExtractUmi(self, log, input_file, out_dir):
        output_file = os.path.join(out_dir, os.path.splitext(os.path.split(input_file)[1])[0] + "_umi.fastq")
        BinaryRunner(BinaryPaths().umi_to_fastq, input_file, output_file).Run(log)
        return output_file

    def CompressFastq(self, log, input_file, out_dir):
        output_file = os.path.join(out_dir, os.path.splitext(os.path.split(input_file)[1])[0] + "_compressed.fastq")
        BinaryRunner(BinaryPaths().ig_trie_compressor, input_file, output_file).Run(log)
        return output_file

    def ConstructGraph(self, log, input_file, out_dir):
        output_file = os.path.join(out_dir, os.path.splitext(os.path.split(input_file)[1])[0] + ".graph")
        BinaryRunner(BinaryPaths().ig_swgraph_construct, input_file, output_file, "-k 6 --tau 1 -A").Run(log)
        return output_file

    def PrintStats(self, log, reads_file, graph_file, output_file):
        BinaryRunner(BinaryPaths().print_graph_decomposition_stats, "", output_file, "-r %s -g %s" % (reads_file, graph_file)).Run(log)

def main():
    log = CreateLogger()
    params = ParseCommandLineParams(log)

    WorkflowRunner().Run(log, params)

if __name__ == '__main__':
    main()
