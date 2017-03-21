import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def rename_headers(germline_file):
    print "Replacing headers in " + germline_file + " file"
    input_records = list(SeqIO.parse(open(germline_file),'fasta'))
    for record in input_records:
        splits = record.id.split('|')
        if len(splits) == 1:
            print "Header " + record.id + " has unknown format!"
            sys.exit(1)
        record.id = splits[1]
        record.description = ""
    print str(len(input_records)) + " records were processed"
    SeqIO.write(input_records, germline_file, "fasta")
    print "FASTA sequences with replaced headers was written to " + germline_file

def main(argv):
    germline_dir = argv[1]
    inner_dirs = os.listdir(germline_dir)
    for inner_dir in inner_dirs:
        inner_dir = os.path.join(germline_dir, inner_dir)
        if not os.path.isdir(inner_dir):
            continue
        inner_dirs2 = os.listdir(inner_dir)
        for inner_dir2 in inner_dirs2:
            inner_dir2 = os.path.join(inner_dir, inner_dir2)
            print "Processing files in directory " + inner_dir2
            germline_files = os.listdir(inner_dir2)
            print "Directory contains " + str(len(germline_files)) + " files"
            for germline_file in germline_files:
                rename_headers(os.path.join(inner_dir2, germline_file))
            print "=================="

if __name__ == '__main__':
    main(sys.argv)