#!/usr/bin/env python2

import os
import sys
from Bio import SeqIO

from argparse import ArgumentParser
from simulate import *


def get_id(id):
    import re
    m = re.match(r"^cluster___(\d+)___size___\d+$", id)

    if m:
        g = m.groups()
        return g[0].strip()
    else:
        return None

assert get_id("") is None

assert get_id("cluster___27833___size___137") == "27833"


if __name__ == "__main__":
    parser = ArgumentParser(description="Convert pRESTO output repertoire to quast format")
    parser.add_argument("-i",
                        dest="input",
                        type=str,
                        help="input FASTA/FASTQ file")
    parser.add_argument("-o",
                        dest="output",
                        type=str,
                        help="output FASTA/FASTQ file")
    parser.add_argument("-r",
                        dest="input_rcm",
                        type=str,
                        help="rcm file for the input FASTA/FASTQ file",
                        default=None)
    parser.add_argument("-R",
                        dest="output_rcm",
                        type=str,
                        help="rcm file for the output FASTA/FASTQ file",
                        default=None)

    args = parser.parse_args()
    print "Command line: %s" % " ".join(sys.argv)

    input_format = idFormatByFileName(args.input)
    output_format = idFormatByFileName(args.output)

    to_drop_ids = set()
    with smart_open(args.input) as fh, smart_open(args.output, "w") as fout:
        for record in SeqIO.parse(fh, input_format):
            if 'N' in record.seq:
                if args.input_rcm is not None:
                    id = get_id(record.description)
                    to_drop_ids.add(id)
            else:
                SeqIO.write(record, fout, output_format)

    if args.input_rcm is not None:
        with open(args.input_rcm) as fin, open(args.output_rcm, "w") as fout:
            for line in fin:
                id_cluster = line.strip().split('\t')
                if len(id_cluster) == 1 or id_cluster[1] not in to_drop_ids:
                    fout.write(line)
                else:
                    fout.write(id_cluster[0] + "\t\n")

    print "Conversion %s -> %s done" % (args.input, args.output)
