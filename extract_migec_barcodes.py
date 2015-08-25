#!/usr/bin/env python

import os
import sys

def LineIsHeader(line):
    return line[0] == '@'

def main():
    if len(sys.argv) != 3:
        print "python extract_migec_barcodes.py assembled_barcodes.fastq output_barcodes.txt"
        sys.exit(1)
    fastq_filename = sys.argv[1]
    if not os.path.exists(fastq_filename):
        print "ERROR: FASTQ file " + fastq_filename + " was not found"
        sys.exit(1)
    in_fhandler = open(fastq_filename, 'r')
    fastq_lines = in_fhandler.readlines()
    barcodes = list()
    for l in fastq_lines:
        if LineIsHeader(l):
            l = l.strip().split(':')
            if len(l) != 3:
                continue
            barcodes.append(l[1])
    print str(len(barcodes)) + " barcodes were extracted from " + fastq_filename
    barcodes_fname = sys.argv[2]
    out_fhandler = open(barcodes_fname, "w")
    for barcode in barcodes:
        out_fhandler.write(barcode + "\n")
    out_fhandler.close()

if __name__ == '__main__':
    main()
