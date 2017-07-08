import os

from collections import OrderedDict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from config.config import config, read_config
from config.parse_input_args import parse_args
from special_utils.os_utils import smart_makedirs

from spots.spots import coldspots, hotspots

from chains.chains import Chains
from mutation_strategies.mutation_strategies import MutationStrategies

from kmer_utilities.filtering_kmers_utilities import filter_by_coverage
from kmer_utilities.kmer_utilities import kmer_index, \
                                          mutative_bases, \
                                          reverse_mut_kmer
from sample_reader.standard_samples import concatenate_kmer_freq_matrices
from mutability_diversity.calculate_mutability_diversity \
    import calculate_mutability_full, \
           calculate_substitution

from fr_cdr_comparison import ttest_pv, mannwhitneyu_pv


def compare_forward_backward_mutations(matrices, figures_dir):
    def compare_forward_backward_mutations_cs(matrices,
                                              chain, strategy,
                                              figures_dir, figure_file,
                                              stat_crit=mannwhitneyu_pv,
                                              sign_level=0.01):
        def threshold_good_fr_cdr(x, axis):
            assert axis == 1
            mut = np.nansum(x[:, [0, 2], :], axis=axis)
            nmut = np.nansum(x[:, [1, 3], :], axis=axis)
            mut, nmut  = mut[:, :, np.newaxis], nmut[:, :, np.newaxis]
            c = np.concatenate((mut, nmut), axis=2)
            return np.min(c, axis=2)

        cov_thresh = 20
        matrices = filter_by_coverage(matrices,
                                      coverage_threshold=cov_thresh,
                                      threshold_function=threshold_good_fr_cdr,
                                      mean_function=np.nanmean)
        muts = calculate_mutability_full(matrices)
        substs = calculate_substitution(matrices)

        muts = muts[:, np.newaxis, :]
        orient_probs = muts * substs

        pvalues = set()
        mutation_ratios = set()
        for i, kmer in enumerate(matrices.kmer_names):
            mut_bases = mutative_bases(kmer)
            k = len(kmer)
            cent_b = kmer[k // 2]
            for j, base in enumerate(mut_bases):
                mut_kmer = reverse_mut_kmer(kmer, base)
                try:
                    back_i = matrices.kmer_names.index(mut_kmer)
                except ValueError:
                    continue
                back_j = mutative_bases(mut_kmer).index(cent_b)
                f_m = np.nanmean(matrices.matrices[i, j + 4, :])
                b_m = np.nanmean(matrices.matrices[back_i, back_j + 4, :])
                if f_m < cov_thresh and b_m < cov_thresh:
                    continue

                forward_probs = orient_probs[i, j, :]
                backward_probs = orient_probs[back_i, back_j]
                pvalues.add(stat_crit(forward_probs, backward_probs,
                                      "two-sided"))
                mutation_ratio = \
                    np.nanmean(forward_probs) / np.nanmean(backward_probs)
                if mutation_ratio < 1.:
                    mutation_ratio = 1. / mutation_ratio
                if not (np.isnan(mutation_ratio) or np.isinf(mutation_ratio)):
                    mutation_ratios.add(mutation_ratio)

        pvalues = np.array(list(pvalues))
        pvalues = np.sort(pvalues)
        pvalues *= np.arange(len(pvalues), 0, -1)
        good_pv = pvalues < sign_level

        mutation_ratios = np.array(list(mutation_ratios))
        if len(mutation_ratios):
            g = sns.distplot(mutation_ratios)
            plt.xlim(xmin=0, xmax=50)
            fig = g.get_figure()
            fig.set_size_inches(4, 3)
            fig.savefig(os.path.join(figures_dir, figure_file),
                        format='pdf', dpi=150)
            plt.close()

        return OrderedDict([
                   ("# well-covered substitutions out of 3072", len(good_pv)),
                   ("% occ with significant difference in forward and backward", np.mean(good_pv)),
                   ("# occ with significant difference in forward and backward", np.sum(good_pv))
               ])


    res = OrderedDict()
    for strategy in MutationStrategies:
        for chain in Chains:
            matr = matrices[strategy][chain]
            figure_file = "distplot_%s_%s.pdf" % (strategy.name, chain.name)
            res[(strategy, chain)] = \
                compare_forward_backward_mutations_cs(
                    matr,
                    strategy, chain,
                    figures_dir, figure_file)

    index = res[(MutationStrategies.NoKNeighbours, Chains.IGH)].keys()
    return pd.DataFrame(res, index=index)
