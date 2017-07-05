#!/usr/bin/env python2

from Bio import SeqIO
import sys


def no_Ns(read):
    return "N" not in read.seq


if __name__ == "__main__":
    input = sys.argv[1]
    if (len(sys.argv) == 3):
        output = sys.argv[2]
    else:
        output = input

    print "Input file: %s" % input
    print "Output file: %s" % output

    with open(input, "rU") as fh:
        reads = list(SeqIO.parse(fh, "fasta"))

    print "Input reads: %d" % len(reads)

    filtered_reads = [read for read in reads if no_Ns(read)]

    print "Output reads: %d" % len(filtered_reads)

    with open(output, "w") as fh:
        SeqIO.write(filtered_reads, fh, "fasta")
