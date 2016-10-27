import os
import sys
from collections import defaultdict

current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = os.path.join(current_dir, os.pardir, os.pardir, os.pardir)
sys.path.append(igrec_dir + "/py/")

from ig_compress_equal_clusters import parse_cluster_mult
from rcm_utils import read_rcm_list


def fix_migec_mixcr_cluster_sizes(input_file, rcm_file, output_file):
    from Bio import SeqIO
    reads = list(SeqIO.parse(input_file, "fasta"))
    rcm = read_rcm_list(rcm_file)
    cluster_size = defaultdict(int)
    cluster_mig_size = defaultdict(int)
    for read_id, cluster in rcm:
        size = int(read_id.split(":")[-1])
        cluster_size[cluster] += size
        cluster_mig_size[cluster] += 1

    with open(output_file, "w") as of:
        for read in reads:
            cluster_id, mult = parse_cluster_mult(read.id)
            assert cluster_mig_size[cluster_id] == mult
            read.id = read.description = "cluster___%s___size___%d" % (cluster_id, cluster_size[cluster_id])
            SeqIO.write(read, of, "fasta")
