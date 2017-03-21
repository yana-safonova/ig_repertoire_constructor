import random
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os

igrec_dir = os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))))

def smart_open(filename, mode="r"):
    import gzip
    import re

    if "w" in mode:
        MODE = "w"
    elif "a" in mode:
        MODE = "a"
    else:
        MODE = "r"

    if re.match(r"^.*\.gz$", filename):
        assert(MODE != "a")
        fh = gzip.open(filename, mode=MODE)
    else:
        fh = open(filename, mode=mode)
    return fh


def get_consensus_dists(all_dists):
    print "reading consensus dists"
    all_dist_file = smart_open(all_dists)
    alld = defaultdict(int)
    d = 0
    # total_dists = 0
    for line in all_dist_file:
        alld[d] = int(line)
        d += 1
        # total_dists += 1
    return alld


def get_dists_inside_clusters(cluster_dists):
    print "reading dists inside clusters"
    cluster_dist_file = smart_open(cluster_dists)
    clusterd = defaultdict(int)
    cluster_n = 0
    while True:
        line = cluster_dist_file.readline()
        if not line:
            break
        size = int(line)
        distn = size * (size - 1) / 2
        dists = [int(cluster_dist_file.readline()) for i in xrange(distn)]
        for d in dists:
            clusterd[d] += 1.0 / distn
        cluster_n += 1

        # if random.randrange(1000) == 5:
        #     break
    return cluster_n, clusterd


def main():
    all_dists = sys.argv[1]
    cluster_dists = sys.argv[2]
    # total_reads = int(sys.argv[3])
    save_file = sys.argv[3]

    alld = get_consensus_dists(all_dists)
    # total_dists = total_reads * (total_reads - 1) / 2
    total_dists = sum(alld.values())
    m = max(alld)
    print "max is", m

    cluster_n, clusterd = get_dists_inside_clusters(cluster_dists)

    print "drawing histograms"
    bin_width = 1
    all_distr = [float(sum([alld[i] for i in alld if i <= d])) / total_dists for d in range(m)]
    print "all distribution", zip(range(len(all_distr)), all_distr)
    dp = sns.distplot(range(m), m / bin_width, hist_kws={'weights': all_distr}, color='blue', kde=False)
    cluster_distr = [float(sum([clusterd[i] for i in clusterd if i >= d])) / cluster_n for d in range(m)]
    print "cluster distribution", zip(range(len(cluster_distr)), cluster_distr)
    dp = sns.distplot(range(m), m / bin_width, hist_kws={'weights': cluster_distr}, color='red', kde=False)
    plt.xlabel("Distance")
    plt.ylabel("% of read pairs")
    plt.title("Distributions of distances between all reads and between reads sharing barcode")
    plt.xlim(0, 100)
    plt.ylim(0, 1)

    # fig = dp.get_figure()
    # fig.savefig(save_file, format=os.path.splitext(save_file)[1][1:])

    from matplotlib.backends.backend_pdf import PdfPages

    pp = PdfPages(save_file)
    pp.savefig()
    pp.close()


if __name__ == "__main__":
    main()