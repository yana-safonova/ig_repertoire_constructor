import sys

from Bio import SeqIO

reads_path = sys.argv[1]
min_length = 100000000
max_length = 0
with open(reads_path) as reads_file:
    records = SeqIO.parse(reads_file, "fasta")
    for record in records:
        length = len(record.seq)
        min_length = min(min_length, length)
        max_length = max(max_length, length)
print "min length: %d, max length: %d" % (min_length, max_length)