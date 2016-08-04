import sys
from collections import defaultdict
from Bio import SeqIO
from ig_compress_equal_clusters import parse_cluster_mult
from rcm_utils import read_rcm_list, rcm2rcmint


def merge():
    fasta_file = sys.argv[1]
    rcm_file = sys.argv[2]
    final_fasta_file = sys.argv[3]
    print "reading rcm file from %s" % rcm_file
    rcm = read_rcm_list(rcm_file)
    cluster_size = defaultdict(int)
    for read_id, cluster_number in rcm:
        cluster_size[cluster_number] += 1
    with open(fasta_file, "r") as in_file, open(final_fasta_file, "w") as out_file:
        for record in SeqIO.parse(in_file, "fasta"):
            cluster_id, mult = parse_cluster_mult(record.description)
            assert cluster_size[cluster_id] > 0
            record.id = record.description = "cluster___%s___size___%d" % (cluster_id, cluster_size[cluster_id])
            SeqIO.write(record, out_file, "fasta")


if __name__ == "__main__":
    merge()