import sys
from collections import defaultdict

from Bio import SeqIO

from rcm_utils import read_rcm_list, rcm2rcmint


def merge():
    fasta_file = sys.argv[1]
    rcm_file = sys.argv[2]
    final_fasta_file = sys.argv[3]
    print "reading rcm file from %s" % rcm_file
    rcm = rcm2rcmint(read_rcm_list(rcm_file))
    cluster_size = defaultdict(int)
    for read_id, cluster_number in rcm:
        cluster_size[cluster_number] += 1
    print "reading repertoire from %s and writing fixed records to %s" % (fasta_file, final_fasta_file)
    with open(fasta_file, "r") as in_file, open(final_fasta_file, "w") as out_file:
        for i, record in enumerate(SeqIO.parse(in_file, "fasta")):
            record.id = record.description = "cluster___%d___size___%d" % (i, cluster_size[i])
            SeqIO.write(record, out_file, "fasta")


if __name__ == "__main__":
    merge()