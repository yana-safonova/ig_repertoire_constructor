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

run_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'

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
    def __init__(self, input_file, output_file, params):
        self.input_file = input_file
        self.output_file = output_file
        self.params = params

class WorkflowRunner:
    def Run(self, log, params):
        umi_fastq = self.ExtractUmi(log, params.input_file, params.tmp_dir)
        umi_compressed = self.CompressFastq(log, umi_fastq, params.tmp_dir)
        umi_graph = self.ConstructGraph(log, umi_compressed, params.tmp_dir)
        self.PrintStats(log, umi_graph, params.stats_file)

    def ExtractUmi(self, log, input_file, out_dir):


def main():
    log = CreateLogger()
    params = ParseCommandLineParams(log)

    WorkflowRunner().Run(log, params)

if __name__ == '__main__':
    main()
