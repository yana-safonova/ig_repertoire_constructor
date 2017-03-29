#!/usr/bin/env python2

from argparse import ArgumentParser
from Bio import SeqIO
import sys


import os
current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = current_dir + "/../"
sys.path.append(igrec_dir + "/py/utils")
sys.path.append(igrec_dir + "/py/pipeline/")
import support
sys.path.append(igrec_dir + "/py/")
from ash_python_utils import idFormatByFileName, smart_open

# from ig_compress_equal_clusters import parse_cluster_mult
def parse_size(s):
    import re

    m = re.match(r"^.*___size___(\d+)$", s)

    if m:
        g = m.groups()
        return int(g[0])
    else:
        return None

assert(parse_size("dsdsfsd___size___10") == 10)


if __name__ == "__main__":
    parser = ArgumentParser(description="Convert repertoire for UCHIME")
    parser.add_argument("input",
                        type=str,
                        help="input FASTA/FASTQ file with abundances in ids")
    parser.add_argument("output",
                        type=str,
                        help="output FASTA/FASTQ file")

    args = parser.parse_args()

    with smart_open(args.input, "r") as fin, smart_open(args.output, "w") as fout:
        for record in SeqIO.parse(fin, idFormatByFileName(args.input)):
            id = str(record.description)
            size = parse_size(id)
            record.id = record.description = id + ";size=%d;" % size
            SeqIO.write(record, fout, idFormatByFileName(args.output))
