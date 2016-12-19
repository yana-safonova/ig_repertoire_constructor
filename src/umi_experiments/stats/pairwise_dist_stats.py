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


def main():
    # rep_file = sys.argv[1]
    # all_dists = sys.argv[2]
    # threads = int(sys.argv[3])
    # max_dist = int(sys.argv[4])
    # print "running dist calculation"
    # os.system("%s %s %s %d %d" % (os.path.join(igrec_dir, "build/release/bin/pairwise_dist_stats"), rep_file, all_dists, threads, max_dist))
    all_dists = sys.argv[1]
    cluster_dists = sys.argv[2]
    save_file = sys.argv[3]

    all_dist_file = smart_open(all_dists)
    # alld = [int(all_dist_file.readline()) for i in xrange(10000)]
    alld = defaultdict(int)
    total_dists = 0
    for line in all_dist_file:
        alld[int(line)] += 1
        total_dists += 1

        # if total_dists > 10000:
        #     break
    # alld = [int(line) for line in all_dist_file]
    print "read %d dists" % len(alld)
    m = max(alld)
    print "max is", m
    # alld = [x for x in alld if x != m]
    print "now only %d dists are left" % len(alld)

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
    print "cluster dists read"

    bin_width = 1
    dp = sns.distplot(range(m), m / bin_width, hist_kws={'weights':[float(alld[d]) / total_dists for d in range(m)]}, color='blue', kde=False)
    dp = sns.distplot(range(m), m / bin_width, hist_kws={'weights':[float(clusterd[d]) / cluster_n for d in range(m)]}, color='red', kde=False)
    plt.xlabel("Distance")
    plt.ylabel("% of read pairs")
    plt.title("Distributions of distances between all reads and between reads sharing barcode")

    # fig = dp.get_figure()
    # fig.savefig(save_file, format=os.path.splitext(save_file)[1][1:])

    from matplotlib.backends.backend_pdf import PdfPages

    pp = PdfPages(save_file)
    pp.savefig()
    pp.close()

if __name__ == "__main__":
    main()