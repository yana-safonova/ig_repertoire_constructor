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
    # intermed_file = sys.argv[2]
    # threads = int(sys.argv[3])
    # max_dist = int(sys.argv[4])
    # print "running dist calculation"
    # os.system("%s %s %s %d %d" % (os.path.join(igrec_dir, "build/release/bin/pairwise_dist_stats"), rep_file, intermed_file, threads, max_dist))
    intermed_file = sys.argv[1]
    save_file = sys.argv[2]

    dist_file = smart_open(intermed_file)
    d = [int(dist_file.readline()) for i in xrange(1000000)]
    # d = [int(line) for line in dist_file]
    print "read %d dists" % len(d)
    m = max(d)
    print "max is", m
    d = [x for x in d if x != m]
    print "now only %d dists are left" % len(d)
    dp = sns.distplot(d, m / 3)
    # plt.show()
    fig = dp.get_figure()
    fig.savefig(save_file, format=os.path.splitext(save_file)[1][1:])

if __name__ == "__main__":
    main()