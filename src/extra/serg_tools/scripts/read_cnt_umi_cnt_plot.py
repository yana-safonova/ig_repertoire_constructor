import os
import sys
from Bio import SeqIO
current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = os.path.join(current_dir, os.pardir, os.pardir, os.pardir, os.pardir)
sys.path.append(igrec_dir)
sys.path.append(igrec_dir + "/py/")
from ig_compress_equal_clusters import parse_cluster_mult
sys.path.append(igrec_dir + "/py/")
from ash_python_utils import smart_open
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def main():
    with smart_open(sys.argv[1], "r") as handle:
        cluster_to_read_cnt = {}
        for record in SeqIO.parse(handle, "fasta"):
            cluster, mult = parse_cluster_mult(record.id)
            cluster_to_read_cnt[cluster] = mult
    with smart_open(sys.argv[2], "r") as handle:
        cluster_to_umi_cnt = {}
        for record in SeqIO.parse(handle, "fasta"):
            cluster, mult = parse_cluster_mult(record.id)
            cluster_to_umi_cnt[cluster] = mult

    read_cnt = np.array([cluster_to_read_cnt[cluster] for cluster in cluster_to_read_cnt.keys()])
    umi_cnt = np.array([cluster_to_umi_cnt[cluster] for cluster in cluster_to_read_cnt.keys()])
    plot = sns.regplot(umi_cnt, read_cnt, fit_reg = False)
    plot.set_ylabel("Read count")
    plot.set_xlabel("Barcode count")
    margin_coef = 0.01
    # ysize = max(read_cnt) - min(read_cnt)
    # ymargin = margin_coef * ysize
    xsize = max(umi_cnt) - min(umi_cnt)
    xmargin = margin_coef * xsize
    ymargin = 0.1
    plt.ylim(1.0 / (1 + ymargin), max(read_cnt) * (1 + ymargin))
    plt.xlim(min(umi_cnt) - xmargin, max(umi_cnt) + xmargin)
    plt.yscale("log", nonposy="clip")
    plt.savefig(os.path.join(os.path.dirname(sys.argv[1]), "read_cnt_to_umi_cnt.png"))
    plt.close()

if __name__ == '__main__':
    main()
