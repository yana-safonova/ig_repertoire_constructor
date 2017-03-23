import os
import sys
from collections import defaultdict

from Bio import SeqIO

current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = os.path.join(current_dir, os.pardir, os.pardir, os.pardir, os.pardir)
sys.path.append(igrec_dir)
sys.path.append(igrec_dir + "/py/")
from igquast_impl import parse_rcm
from ig_compress_equal_clusters import parse_cluster_mult

reference_path = sys.argv[1]
reference_rcm_path = sys.argv[2]
# barigrec_rcm_path = sys.argv[3]

large_cluster_id = "not found"
with open(reference_path) as reference_file:
    records = SeqIO.parse(reference_file, "fasta")
    for record in records:
        id, mult = parse_cluster_mult(record.id)
        if mult > 2000:
            large_cluster_id = id
            break
print "large cluster: ", large_cluster_id

reference_rcm = parse_rcm(reference_rcm_path)
cluster = [key for key in reference_rcm.keys() if reference_rcm[key] == large_cluster_id]
print len(cluster)

rev = defaultdict(int)
for key, value in reference_rcm.iteritems():
    rev[value] += 1

print max(rev.values())