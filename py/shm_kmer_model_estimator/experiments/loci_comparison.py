import os.path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['figure.figsize'] = 15, 10

from collections import OrderedDict

from chains.chains import Chains
from mutation_strategies.mutation_strategies import MutationStrategies

from kmer_utilities.filtering_kmers_utilities import filter_by_coverage
from mutability_diversity.calculate_mutability_diversity \
    import calculate_mutability_full

from scipy.stats.mstats import kruskalwallis


def kruskalwallis_pv(*args):
    return kruskalwallis(*args).pvalue


def compare_loci(matrices, figures_dir):
    def compare_loci_strategy(matrices, strategy,
                              figures_dir, figure_file,
                              stat_crit=kruskalwallis_pv,
                              sign_lev=0.01):
        matr_h = matrices[Chains.IGH]
        matr_l = matrices[Chains.IGL]
        matr_k = matrices[Chains.IGK]

        muts_h = calculate_mutability_full(matr_h)
        muts_l = calculate_mutability_full(matr_l)
        muts_k = calculate_mutability_full(matr_k)

        pvalues = []
        means_h, means_l, means_k = [], [], []
        for i, kmer in enumerate(matr_h.kmer_names):
            mut_h, mut_l, mut_k = muts_h[i, :], muts_l[i, :], muts_k[i, :]

            mut_h = mut_h[~np.isnan(mut_h)]
            mut_k = mut_k[~np.isnan(mut_k)]
            mut_l = mut_l[~np.isnan(mut_l)]
            if len(mut_h) * len(mut_k) * len(mut_l) == 0:
                continue

            if np.all(mut_h == mut_h[0]) or \
                np.all(mut_k == mut_k[0]) or \
                np.all(mut_l == mut_l[0]):
                continue

            pvalues.append(stat_crit(mut_h, mut_l, mut_k))
            means_h.append(np.nanmean(mut_h))
            means_l.append(np.nanmean(mut_l))
            means_k.append(np.nanmean(mut_k))

        pvalues = np.array(list(pvalues))
        pvalues = np.sort(pvalues)
        pvalues *= np.arange(len(pvalues), 0, -1)
        good_pv = pvalues < sign_lev

        means_h = np.array(means_h)
        means_l = np.array(means_l)
        means_k = np.array(means_k)

        means = np.vstack([means_h, means_l, means_k]).T
        means = pd.DataFrame(means,
                             columns=["IGH", "IGL", "IGK"])
        means = pd.melt(means, var_name='loci', value_name='mutability')
        g = sns.violinplot(x="loci", y="mutability", data=means, 
                           inner="quartile")
        plt.ylim([0, 0.55])
        fig = g.get_figure()
        fig.savefig(os.path.join(figures_dir, figure_file),
                    format='pdf', dpi=150)
        plt.close()
        return OrderedDict([
                   ("# well-covered kmers", len(good_pv)),
                   ("% kmers with significant mut_FR < mut_CDR", np.mean(good_pv)),
                   ("# kmers with significant mut_FR < mut_CDR", np.sum(good_pv))
               ])

    results = OrderedDict()
    for strategy in MutationStrategies:
        figure_file = "violinplot_%s.pdf" % strategy.name
        results[strategy] = \
            compare_loci_strategy(matrices[strategy],
                                  strategy,
                                  figures_dir, figure_file)
    index = results[MutationStrategies.NoKNeighbours].keys()
    return pd.DataFrame(results, index=index)
