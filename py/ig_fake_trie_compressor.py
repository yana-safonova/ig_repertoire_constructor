#!/usr/bin/env python2

from argparse import ArgumentParser
from Bio import SeqIO

import sys


import os
current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = current_dir + "/../"
sys.path.append(igrec_dir + "/py/utils")
sys.path.append(igrec_dir + "/py/pipeline/")
sys.path.append(igrec_dir + "/py/")
from ash_python_utils import idFormatByFileName, smart_open

if __name__ == "__main__":
    parser = ArgumentParser(description="Fake ig_trie_compressor")
    parser.add_argument("--input", "-i",
                        required=True,
                        type=str,
                        help="input FASTA/FASTQ file with abundances in ids")
    parser.add_argument("--output", "-o",
                        required=True,
                        type=str,
                        help="output FASTA/FASTQ file")
    parser.add_argument("--id-map", "-m",
                        type=str,
                        default="",
                        help="map file name; empty (default) for non-producing")

    args = parser.parse_args()

    print "Fake trie_compressor started"
    read_count = 0
    with smart_open(args.input, "r") as fin, smart_open(args.output, "w") as fout:
        for record in SeqIO.parse(fin, idFormatByFileName(args.input)):
            id = str(record.description)
            record.id = record.description = id + "__size__1"
            SeqIO.write(record, fout, idFormatByFileName(args.output))
            read_count += 1

    if args.id_map:
        with smart_open(args.id_map, "w") as f_id_map:
            for i in xrange(read_count):
                f_id_map.write("%d\n" % i)

    print "Fake trie_compressor finished"
