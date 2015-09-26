#!/usr/bin/env python2

from argparse import ArgumentParser
from barcode_cleaner import smart_open

if __name__ =="__main__":
    parser = ArgumentParser("Replace spaces by underlines")
    parser.add_argument("input",
                        type=str,
                        help="input FASTQ file")
    parser.add_argument("output",
                        type=str,
                        help="output FASTQ file")


    args = parser.parse_args()

    with smart_open(args.input, "r") as fi:
        with smart_open(args.output, "w") as fo:
            for l in fi:
                l = l.replace(" ", "_")
                fo.write(l)
