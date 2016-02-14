#!/usr/bin/env python2

from argparse import ArgumentParser
from Bio import SeqIO
import sys
import os
import os.path
import tempfile


from ig_remove_low_abundance_reads import smart_open, parse_abundance



if __name__ == "__main__":
    parser = ArgumentParser(description="Compress equal clusters")
    parser.add_argument("input",
                        type=str,
                        help="input FASTA file")
    parser.add_argument("output",
                        type=str,
                        help="output FASTA file")
    parser.add_argument("--tmp-file", "-T",
                        type=str,
                        default="",
                        help="temporary file for ig_trie_compressor output (default: EMPTY")
    parser.add_argument("--map-file", "-m",
                        type=str,
                        default="",
                        help="map file (default: EMPTY")

    args = parser.parse_args()

    if (not args.tmp_file):
        args.tmp_file = tempfile.mkstemp(suffix=".fa", prefix="igrec_")[1]
    if (not args.map_file):
        args.map_file = tempfile.mkstemp(suffix=".map", prefix="igrec_")[1]

    print "Command line: %s" % " ".join(sys.argv)

    home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/../../'
    run_trie_compressor = os.path.join(home_directory, 'build/release/bin/./ig_trie_compressor')
    cmd_line = "%s -i %s -o %s -m %s" % (run_trie_compressor, args.input, args.tmp_file, args.map_file)
    os.system(cmd_line)

    with open(args.map_file) as fin:
        targets = [int(line.strip()) for line in fin]

    with open(args.input) as fin:
        mults = [parse_abundance(record.description) for record in SeqIO.parse(fin, "fasta")]

    final_mults = [0] * (max(targets) + 1) if len(targets) else []
    for target, mult in zip(targets, mults):
        final_mults[target] += mult

    with smart_open(args.tmp_file, "r") as fin, smart_open(args.output, "w") as fout:
        for i, record in enumerate(SeqIO.parse(fin, "fasta")):
            record.id = record.description = "cluster___%d___size___%d" % (i + 1, final_mults[i])
            SeqIO.write(record, fout, "fasta")

    print "Remove temporary files: %s %s" % (args.tmp_file, args.map_file)
    os.remove(args.tmp_file)
    os.remove(args.map_file)
