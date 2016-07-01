#!/usr/bin/env python2

from argparse import ArgumentParser
from Bio import SeqIO
import sys
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams.update({'font.size': 16})



igrc_path = "/home/ashlemov/Git/ig_repertoire_constructor"

trie_comp = "%s/build/release/bin/ig_trie_compressor " % igrc_path
ig_matcher_bin = "%s/build/release/bin/ig_matcher " % igrc_path

def parse_multiplicity(s):
    import re

    m = re.match(r".*_multiplicity_(\d+)_", s)

    if m:
        g = m.groups()
        return int(g[0])
    else:
        return 1

def parse_size(s):
    import re

    m = re.match(r".*___size___(\d+)", s)

    if m:
        g = m.groups()
        return int(g[0])
    else:
        return 1

def mult2mult(clustering_fa, reference_fa, match_filename, simulator=False):
    with open(clustering_fa) as f:
        clustering_mults = [parse_size(rec.id) for rec in SeqIO.parse(f, "fasta")]

    with open(reference_fa) as f:
        if simulator:
            reference_mults = [parse_multiplicity(rec.id) for rec in SeqIO.parse(f, "fasta")]
        else:
            reference_mults = [parse_size(rec.id) for rec in SeqIO.parse(f, "fasta")]

    unmatched_references = set(range(len(reference_mults)))
    result = []

    with open(match_filename) as f:
        for i, line in enumerate(f):
            line = line.strip()
            cl_mult = clustering_mults[i]
            ref_mult = 0
            if line:
                dist = -int(line.split()[0])
                if dist == 0:
                    neibs = map(int, line.split()[1:])
                    if len(neibs):
                        ref_mult = sum(reference_mults[j] for j in neibs)
                        unmatched_references.difference_update(neibs)

            result.append((cl_mult, ref_mult))

    for i in unmatched_references:
        result.append((0, reference_mults[i]))

    return result

def match_fas(clustering_fa, reference_fa, ig_matcher_bin=ig_matcher_bin,
              simulator=False):
    from os import system

    # TODO Use proper tmps
    system("%s -i %s -I %s -o m1.match -O m2.match -k 10 --tau 1" % (ig_matcher_bin,
                                                                     reference_fa,
                                                                     clustering_fa))

    m2m = mult2mult(clustering_fa, reference_fa, "m2.match",
                    simulator=simulator)

    clustering_mults = [clust_mult for clust_mult, reference_mults in m2m]
    reference_mults = [reference_mults for clust_mult, reference_mults in m2m]

    reference_mults_sorted = sorted(reference_mults)
    def ref_count(limit):
        return how_many_greater_or_equal(limit, reference_mults_sorted)

    clustering_mults_sorted = sorted(clustering_mults)
    def cluster_count(limit):
        return how_many_greater_or_equal(limit, clustering_mults_sorted)

    clustering_paired_mults, reference_paired_mults = zip(*m2m)
    minimum_paired_mults = [min(clust_mult, ref_mult) for clust_mult, ref_mult in m2m]

    minimum_paired_mults_sorted = sorted(minimum_paired_mults)
    def sensitivity(limit):
        assert limit > 0
        d = how_many_greater_or_equal(limit, minimum_paired_mults_sorted)
        dd = ref_count(limit)
        if dd == 0:
            return 0.
        else:
            return float(d) / float(dd)

    def specificity(limit):
        assert limit > 0
        d = how_many_greater_or_equal(limit, minimum_paired_mults_sorted)
        dd = cluster_count(limit)
        if dd == 0:
            return 0.
        else:
            return float(d) / float(dd)

    def fdr(limit):
        return 1.0 - specificity(limit)


    def median_rate(limit):
        # TODO Speedup it
        import numpy as np
        assert limit > 0
        rates = [float(clust_mult) / float(ref_mult) for clust_mult, ref_mult in m2m if min(ref_mult, clust_mult) >= limit]

        return np.median(rates)

    class Res:
        pass

    res = Res()
    res.sensitivity, res.specificity, res.fdr = sensitivity, specificity, fdr

    res.median_rate = median_rate
    res.m2m = m2m

    return res

def how_many_greater_or_equal(limit, x):
    import bisect
    return len(x) - bisect.bisect_left(x, limit)

assert how_many_greater_or_equal(4, [4, 4, 4, 4]) == 4
assert how_many_greater_or_equal(4, [4]) == 1
assert how_many_greater_or_equal(4, []) == 0
assert how_many_greater_or_equal(4, [1, 2, 3, 4, 5]) == 2


res = match_fas("/home/ashlemov/tmp/age3ideal/out/final_repertoire.fa", "/home/ashlemov/tmp/age3ideal/repertoire_compressed.fa")

def plot_mult2mult(res,
                   out="mult2mult",
                   uplimit=1600,
                   bottomlimit=5,
                   title=""):
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import seaborn as sns
    matplotlib.rcParams.update({'font.size': 16})
    sns.set_style("darkgrid")

    f, ax = plt.subplots(figsize=(6, 6))

    m2m = res.m2m
    m2m = [_ for _ in m2m if min(*_) >= bottomlimit]

    uplimit = 1600
    m2m = [_ for _ in m2m if max(*_) < uplimit]


    cluster_mults, reference_mults = zip(*m2m)
    g = sns.JointGrid(x=np.array(reference_mults), y=np.array(cluster_mults),
                    xlim=(0, uplimit), ylim=(0, uplimit),
                    ratio=5).set_axis_labels("Reference cluster size", "Estimated cluster size")
    # g = g.plot_joint(sns.plt.scatter, color="m")
    g = g.plot_joint(sns.plt.scatter)
    g.plot_marginals(sns.distplot, kde=False, color=".5")

    ax = g.ax_joint
    ax.plot([0, 10**10], [0, 10**10], "--")
    ax.plot([0, 10**10], [0, 10**10 * res.median_rate(5)])

    if title:
        plt.title(title)

    plt.savefig(out + ".png")

    pp = PdfPages(out + ".pdf")
    pp.savefig()
    pp.close()



plot_mult2mult(res)



def plot_fdr_sensitivity(res,
                         out="fdr_sensitivity",
                         title=""):
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import seaborn as sns
    matplotlib.rcParams.update({'font.size': 16})
    sns.set_style("darkgrid")

    f, ax = plt.subplots(figsize=(6, 6))

    limits = range(1, 51)

    fdrs = np.array(map(res.fdr, limits))
    sensitivities = np.array(map(res.sensitivity, limits))

    sns.set_style("darkgrid")
    plt.plot(sensitivities, fdrs, "b-")
    plt.xlim((0, 1))
    plt.ylim((0, 1))
    mask = np.array([0, 2, 4, 49])
    plt.plot(sensitivities[mask], fdrs[mask], "bo")
    plt.ylabel("FDR")
    plt.xlabel("Sensitivity")

    plt.annotate('limit = 1', xy=(sensitivities[0], fdrs[0]),
                arrowprops=dict(facecolor='black', shrink=0.05),
                xytext=(sensitivities[0] - 0.25, fdrs[0] + 0.05))

    plt.annotate('limit = 3', xy=(sensitivities[2], fdrs[2]),
                arrowprops=dict(facecolor='black', shrink=0.05),
                xytext=(sensitivities[2] + 0.05, fdrs[2]-0.075))

    plt.annotate('limit = 5', xy=(sensitivities[4], fdrs[4]),
                arrowprops=dict(facecolor='black', shrink=0.05),
                xytext=(sensitivities[4] - 0.15, fdrs[4] + 0.15))

    plt.annotate('limit = 50', xy=(sensitivities[49], fdrs[49]),
                arrowprops=dict(facecolor='black', shrink=0.05),
                xytext=(sensitivities[49] + 0.1, fdrs[49] + 0.1))

    if title:
        plt.title(title)

    plt.savefig(out + ".png")

    pp = PdfPages(out + ".pdf")
    pp.savefig()
    pp.close()


plot_fdr_sensitivity(res)

#
# from rand import make_ideal_rcm, rcm_vs_rcm
#
# make_ideal_rcm("final_repertoire.rcm", "ideal.rcm")
#
# print rcm_vs_rcm("final_repertoire.rcm", "ideal.rcm")

from rand import error_profile
def error_prof(dirname, **kwargs):
    from rand import error_profile
    return error_profile(dirname + "/final_repertoire.rcm",
                         dirname + "/vj_finder/cleaned_reads.fa",
                         dirname + "/final_repertoire.fa",
                         **kwargs)



# ep = error_prof("/home/ashlemov/tmp/age3test/out_new_fees/", limit=10)
# ep = error_prof("/home/ashlemov/tmp/age3ideal/out/")
# ep = error_profile("/home/ashlemov/tmp/age3ideal/out/final_repertoire.rcm",
#                    "/home/ashlemov/tmp/age3ideal/vj_finder/age_3_good.fq",
#                    "/home/ashlemov/tmp/age3ideal/out/final_repertoire.fa")
# ep_ideal = error_profile("/home/ashlemov/tmp/age3ideal/barcodes_3.rcm",
#                          "/home/ashlemov/tmp/age3ideal/out7/vj_finder/cleaned_reads.fa",
#                          "/home/ashlemov/tmp/age3ideal/repertoire_3.fa")
ep_ideal = error_profile("/home/ashlemov/tmp/age3ideal/barcodes_3.rcm",
                         "/home/ashlemov/tmp/age3ideal/vj_finder/age_3_good.fq",
                         "/home/ashlemov/tmp/age3ideal/repertoire_3.fa")

ep_over7 = error_profile("/home/ashlemov/tmp/age3ideal/out7/final_repertoire.rcm",
                         # "/home/ashlemov/tmp/age3ideal/vj_finder/age_3_good.fq",
                         "/home/ashlemov/tmp/age3ideal/out7/vj_finder/cleaned_reads.fa",
                         "/home/ashlemov/tmp/age3ideal/out7/final_repertoire.fa")

ep_usual = error_profile("/home/ashlemov/tmp/age3ideal/out/final_repertoire.rcm",
                         # "/home/ashlemov/tmp/age3ideal/vj_finder/age_3_good.fq",
                         "/home/ashlemov/tmp/age3ideal/out/vj_finder/cleaned_reads.fa",
                         "/home/ashlemov/tmp/age3ideal/out/final_repertoire.fa")

print ep_ideal.error_rate
print ep_usual.error_rate
print ep_over7.error_rate

def plot_error_profile(ep,
                       out="error_profile",
                       title="",
                       annotate_cdrs=False):
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import seaborn as sns
    matplotlib.rcParams.update({'font.size': 16})
    sns.set_style("darkgrid")

    f, ax = plt.subplots(figsize=(6, 6))
    sns.distplot(ep.errors01).set(xlim=(0, 1))

    if annotate_cdrs:
        plt.annotate('CDR1', xy=(0.15, 2.25))
        plt.annotate('CDR2', xy=(0.39, 1.75))
        plt.annotate('CDR3', xy=(0.85, 0.8))

    if title:
        plt.title(title)

    plt.savefig(out + ".png")

    pp = PdfPages(out + ".pdf")
    pp.savefig()
    pp.close()




plot_error_profile(ep_ideal, out="good_ep")
plot_error_profile(ep_over7, out="overcorrected", annotate_cdrs=True)



def plot_errors_in_reads_dist(errors_in_read,
                              out="ErrorsInReadDistribution",
                              title="",
                              max_val = 5):
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import seaborn as sns
    matplotlib.rcParams.update({'font.size': 16})
    sns.set_style("darkgrid")

    lam = np.mean(errors_in_read)

    np.random.seed(0)
    f, ax = plt.subplots(figsize=(6, 6))


    bins = np.array(range(max_val + 1)) - 0.5
    a_heights, a_bins = np.histogram(errors_in_read, bins=bins)
    x = np.random.poisson(lam, len(errors_in_read))
    b_heights, b_bins = np.histogram(x, bins=a_bins)

    width = (a_bins[1] - a_bins[0])/3

    actual = ax.bar(a_bins[:-1]+width/2, a_heights, width=width, facecolor='cornflowerblue',
                    label="Actual frequencies")

    pois_model_name = "Poisson distribution ($\hat{\lambda} = %.2f$)" % lam
    poiss = ax.bar(b_bins[:-1]+width+width/2, b_heights, width=width, facecolor='seagreen',
                label=pois_model_name)
    plt.legend(handles=[actual, poiss])
    plt.xlim(-width, max_val + width)

    if title:
        plt.title(title)

    plt.savefig(out + ".png")

    pp = PdfPages(out + ".pdf")
    pp.savefig()
    pp.close()

plot_errors_in_reads_dist(ep_usual.errors_in_read)
