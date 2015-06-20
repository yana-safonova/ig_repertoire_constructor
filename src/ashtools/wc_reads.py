#!/usr/bin/env python2

from Bio import SeqIO
import argparse


def ilen(iterable):
    return sum(1 for i in iterable)


def determine_format_by_ext(fname):
    import os
    fileName, fileExtension = os.path.splitext(fname)

    if fileExtension in [".fq", ".fastq"]:
        return "fastq"
    elif fileExtension in [".fa", ".fasta"]:
        return "fasta"
    else:
        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count the number of reads in each file (like `wc -l`)")
    parser.add_argument("files", metavar="files", type=str, nargs='+',
                        help="files with reads (FASTA or FASTQ)")

    args = parser.parse_args()
    for fname in args.files:
        with open(fname, "rU") as fh:
            l = ilen(SeqIO.parse(fh, determine_format_by_ext(fname)))
            print("%d\t%s" % (l, fname))
