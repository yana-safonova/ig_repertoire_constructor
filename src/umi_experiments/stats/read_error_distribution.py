import random
from collections import defaultdict
import matplotlib
from os.path import splitext

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


def get_error_cnt_distr(distr_file_name):
    print "reading read error distribution"
    distr_file = smart_open(distr_file_name)
    distr = [int(line) for line in distr_file]
    return distr


def main():
    distr_file_name = sys.argv[1]
    save_file = sys.argv[2]

    distr = get_error_cnt_distr(distr_file_name)

    print "drawing histograms"
    bin_width = 1
    max = len(distr) - 1
    norm_distr = distr[0: max + 1]
    total_reads = sum(norm_distr)
    norm_distr = [float(x) / total_reads for x in norm_distr]
    dp = sns.distplot(range(max + 1), (max + 1) / bin_width, hist_kws={'weights': norm_distr}, color='blue', kde=False)
    plt.xlabel("Error count")
    plt.ylabel("% of reads")
    plt.title("Distribution of amplification errors per read")
    plt.xlim(0, 25)
    # plt.ylim(0, 1)

    # fig = dp.get_figure()
    # fig.savefig(save_file, format=os.path.splitext(save_file)[1][1:])

    format = splitext(save_file)[-1][1:]

    if 'png' in format:
        plt.savefig(save_file, format='png')

    if 'pdf' in format:
        from matplotlib.backends.backend_pdf import PdfPages
        pp = PdfPages(save_file)
        pp.savefig()
        pp.close()

    # from matplotlib.backends.backend_pdf import PdfPages
    #
    # pp = PdfPages(save_file)
    # pp.savefig()
    # pp.close()


if __name__ == "__main__":
    main()