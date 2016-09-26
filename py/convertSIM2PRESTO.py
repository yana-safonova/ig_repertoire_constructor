#!/usr/bin/env python2

import os
import sys
from Bio import SeqIO

from argparse import ArgumentParser
from simulate import *


def convert_id(id, primer="NA"):
    import re
    m = re.match(r"^(.*)UMI:([ACGTN]+):?(.*)$", id)

    if m:
        g = m.groups()
        id = g[0].strip()
        barcode = g[1].strip()
        return id + "|PRIMER=%s|BARCODE=%s" % (primer, barcode)
    else:
        return None

assert convert_id("M01691:10:000000000-A7F7L:1:1101:18193:1582_1:N:0:1_UMI:GCGGAATATTTCA:DF77+,E<FEGGF") == \
        "M01691:10:000000000-A7F7L:1:1101:18193:1582_1:N:0:1_|PRIMER=NA|BARCODE=GCGGAATATTTCA"

assert convert_id("34014_mutated_from_4083_UMI:TTTGTTCAATCTCG") == \
        "34014_mutated_from_4083_|PRIMER=NA|BARCODE=TTTGTTCAATCTCG"


if __name__ == "__main__":
    parser = ArgumentParser(description="Convert AGE dataset to PRESTO")
    parser.add_argument("input",
                        type=str,
                        help="input FASTA/FASTQ file")
    parser.add_argument("output",
                        type=str,
                        help="output FASTA/FASTQ file")
    parser.add_argument("--primer", "-p",
                        type=str,
                        default="NA",
                        help="primer name to add")

    args = parser.parse_args()
    print "Command line: %s" % " ".join(sys.argv)

    input_format = idFormatByFileName(args.input)
    output_format = idFormatByFileName(args.output)

    with smart_open(args.input) as fh, smart_open(args.output, "w") as fout:
        for record in SeqIO.parse(fh, input_format):
            record.id = record.description = convert_id(record.description, args.primer)
            SeqIO.write(record, fout, output_format)

    print "Conversion %s -> %s done" % (args.input, args.output)
