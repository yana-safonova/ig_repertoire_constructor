import os
import sys

from Bio import SeqIO

current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = os.path.join(current_dir, os.pardir, os.pardir, os.pardir)
sys.path.append(igrec_dir + "/py/")

from ig_compress_equal_clusters import parse_cluster_mult

threshold = int(sys.argv[3])

with open(sys.argv[1]) as constructed_file:
    records = SeqIO.parse(constructed_file, "fasta")
    ref_recs = sorted([pair for pair in [(rec.seq, parse_cluster_mult(rec.id)[1]) for rec in records] if pair[1] >= threshold])
with open(sys.argv[2]) as constructed_file:
    records = SeqIO.parse(constructed_file, "fasta")
    const_recs = sorted([(rec.seq, parse_cluster_mult(rec.id)[1]) for rec in records])

i = 0
j = 0
while i < len(ref_recs) and j < len(const_recs):
    if ref_recs[i][0] < const_recs[j][0]:
        i += 1
    elif ref_recs[i][0] > const_recs[j][0]:
        j += 1
    else:
        print float(const_recs[j][1]) / float(ref_recs[i][1]), const_recs[j][1], ref_recs[i][1]
        i += 1
        j += 1

