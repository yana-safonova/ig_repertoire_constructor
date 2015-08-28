import os
import sys
import shutil

from Bio import SeqIO
import sys

def PrepareOutputDir(output_dir):
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

def main():
    if len(sys.argv) != 4:
        print "ERROR: incorrect input parameters"
        print "python create_fastq_barcodes.py reads.fastq barcodes.rcm output_dir"
        sys.exit(1)
    fastq_fname = sys.argv[1]
    result_records = []
    for record in SeqIO.parse(open(fastq_fname), 'fastq'):
        result_records.append(record)
    print str(len(result_records)) + " reads were extracted from " + fastq_fname

    barcodes_fname = sys.argv[2]
    read_barcodes = dict()
    fhandler = open(barcodes_fname, "r")
    for line in fhandler:
        splits = line.strip().split()
        read_barcodes[splits[0]] = splits[1]

    barcode_reads = dict()
    for record in result_records:
        read_name = record.id
        barcode = read_barcodes[read_name]
        if not barcode in barcode_reads:
            barcode_reads[barcode] = list()
        barcode_reads[barcode].append(record)

    output_dir = sys.argv[3]
    PrepareOutputDir(output_dir)
    for barcode in barcode_reads:
        records = barcode_reads[barcode]
        barcode_fname = os.path.join(output_dir, barcode + "_size_" + str(len(record)) + ".txt")
        SeqIO.write(records, open(barcode_fname, 'w'), 'fastq')
        print "Barcode " + barcode + " was written to " + barcode_fname