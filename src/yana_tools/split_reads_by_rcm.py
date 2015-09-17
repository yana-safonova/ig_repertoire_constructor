import os
import sys
import shutil

from Bio import SeqIO
import sys

def PrepareOutputDir(output_dir):
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

def GetReaderDescriptor(reads_fname):
    fa_suffices = ['fa', 'fasta']
    if reads_fname[len(reads_fname) - len(fa_suffices[0]):] == fa_suffices[0] or \
                    reads_fname[len(reads_fname) - len(fa_suffices[1])] == fa_suffices[1]:
        return "fasta"
    fq_suffices = ['fq', 'fastq']
    if reads_fname[len(reads_fname) - len(fq_suffices[0]):] == fq_suffices[0] or \
                    reads_fname[len(reads_fname) - len(fq_suffices[1])] == fq_suffices[1]:
        return "fastq"
    print "ERROR: Reads " + reads_fname + " has wrong description"
    sys.exit(1)

def main():
    if len(sys.argv) != 4:
        print "ERROR: incorrect input parameters"
        print "python create_fastq_rcm.py reads.fast(q/a) read_map.rcm output_dir"
        sys.exit(1)
    reads_fname = sys.argv[1]
    result_records = []
    print "Reads filename: " + reads_fname
    description = GetReaderDescriptor(reads_fname)
    for record in SeqIO.parse(open(reads_fname), description):
        result_records.append(record)
    print str(len(result_records)) + " reads were extracted from " + reads_fname

    rcm_fname = sys.argv[2]
    print "RCM filename: " + rcm_fname
    read_barcodes = dict()
    fhandler = open(rcm_fname, "r")
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
    print "Output directory: " + output_dir
    PrepareOutputDir(output_dir)
    for barcode in barcode_reads:
        records = barcode_reads[barcode]
        if len(records) < 100:
            continue
        barcode_fname = os.path.join(output_dir, barcode + "_size_" + str(len(records)) + "." + description)
        SeqIO.write(records, open(barcode_fname, 'w'), description)
        print "Cluster " + barcode + " was written to " + barcode_fname

if __name__ == '__main__':
    main()
