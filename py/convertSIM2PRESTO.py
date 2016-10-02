#!/usr/bin/env python2

import os
import sys
from Bio import SeqIO

from argparse import ArgumentParser
from simulate import *


def convert_id(id, index):
    import re
    m = re.match(r"^.*DUPCOUNT=(\d+)$", id)

    if m:
        g = m.groups()
        id = "cluster___%d___size___%d" % (index, int(g[0].strip()))
        return id
    else:
        return None

assert convert_id("", 2345) is None

assert convert_id(">1TTTTCAATCGACGGA|PRCONS=NA|PRFREQ=1.0|CONSCOUNT=30|DUPCOUNT=3", 2345) == \
        "cluster___2345___size___3"

if __name__ == "__main__":
    parser = ArgumentParser(description="Convert pRESTO output repertoire to quast format")
    parser.add_argument("input",
                        type=str,
                        help="input FASTA/FASTQ file")
    parser.add_argument("output",
                        type=str,
                        help="output FASTA/FASTQ file")

    args = parser.parse_args()
    print "Command line: %s" % " ".join(sys.argv)

    input_format = idFormatByFileName(args.input)
    output_format = idFormatByFileName(args.output)

    current = 0
    with smart_open(args.input) as fh, smart_open(args.output, "w") as fout:
        for record in SeqIO.parse(fh, input_format):
            record.id = record.description = convert_id(record.description, current)
            SeqIO.write(record, fout, output_format)
            current += 1

    print "Conversion %s -> %s done" % (args.input, args.output)
