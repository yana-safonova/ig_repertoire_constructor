import os.path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from collections import OrderedDict

from chains.chains import Chains
from mutation_strategies.mutation_strategies import MutationStrategies

from kmer_utilities.filtering_kmers_utilities import filter_by_coverage
from mutability_diversity.calculate_mutability_diversity \
    import calculate_mutability_fr, \
           calculate_mutability_cdr

from scipy.stats import ttest_ind, mannwhitneyu


def ttest_pv(a, b):
    return ttest_ind(a, b, equal_var=False).pvalue[0]


def mannwhitneyu_pv(a, b, alternative="less"):
    return mannwhitneyu(a, b, alternative=alternative).pvalue


def compare_fr_cdr(matrices, figures_dir):
    def compare_fr_cdr_by_strategy_chain(matrices,
                                         chain, strategy,
                                         figures_dir, figure_file,
                                         stat_crit=mannwhitneyu_pv,
                                         sign_lev=0.01):
        def threshold_good_fr_cdr(x, axis):
            return np.nanmin(x[:, :4], axis=axis)

        matrices = filter_by_coverage(matrices,
                                      coverage_threshold=20,
                                      threshold_function=threshold_good_fr_cdr,
                                      mean_function=np.nanmean)

        mut_frs = calculate_mutability_fr(matrices)
        mut_cdrs = calculate_mutability_cdr(matrices)

        pvalues = []
        means_fr, means_cdr = [], []
        for i, kmer in enumerate(matrices.kmer_names):
            mut_fr, mut_cdr = mut_frs[i, :], mut_cdrs[i, :]
            pvalues.append(stat_crit(mut_fr, mut_cdr))
            means_fr.append(np.nanmean(mut_fr))
            means_cdr.append(np.nanmean(mut_cdr))

        pvalues = np.array(pvalues)
        pvalues = np.sort(pvalues)
        pvalues *= np.arange(len(pvalues), 0, -1)
        good_pv = pvalues < sign_lev

        means_fr, means_cdr = np.array(means_fr), np.array(means_cdr)
        means = np.vstack([means_fr, means_cdr]).T
        means = pd.DataFrame(means,
                             columns=["FR", "CDR"])
        means = pd.melt(means, var_name='region', value_name='mutability')
        means["all"] = ""
        if not means.empty:
            g = sns.violinplot(x="all", y="mutability", hue="region",
                               data=means, inner="quartile", split=True)
            plt.ylim([0, 0.55])
            fig = g.get_figure()
            fig.set_size_inches(4, 3)
            fig.savefig(os.path.join(figures_dir, figure_file + "pdf"),
                        format='pdf', dpi=150)
            fig.savefig(os.path.join(figures_dir, figure_file + "png"),
                        format='png', dpi=150)
            plt.close()

        return OrderedDict([
                   ("# well-covered kmers", len(matrices.kmer_names)),
                   ("% kmers with significant mut_FR < mut_CDR", np.mean(good_pv)),
                   ("# kmers with significant mut_FR < mut_CDR", np.sum(good_pv))
               ])

    results = OrderedDict()
    for strategy in MutationStrategies:
        for chain in Chains:
            figure_file = "violinplot_%s_%s." % (strategy.name, chain.name)
            results[(strategy, chain)] = \
                compare_fr_cdr_by_strategy_chain(matrices[strategy][chain],
                                                 chain, strategy,
                                                 figures_dir, figure_file)
    index = results[(MutationStrategies.NoKNeighbours, Chains.IGH)].keys()
    return pd.DataFrame(results, index=index)
