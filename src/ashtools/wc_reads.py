#!/usr/bin/env python2

from Bio import SeqIO
import argparse

def ilen(iterable):
    return sum(1 for i in iterable)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count the number of reads in each file (like wc)")
    parser.add_argument("files", metavar="files", type=str, nargs='+',
                        help="files with reads")

    args = parser.parse_args()
    for fname in args.files:
        with open(fname, "rU") as fh:
            l = ilen(SeqIO.parse(fh, "fastq"))
            print("%d\t%s" % (l, fname))
