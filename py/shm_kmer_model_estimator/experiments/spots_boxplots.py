import os

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['figure.figsize'] = 15, 10

from config.config import config, read_config
from config.parse_input_args import parse_args
from special_utils.os_utils import smart_makedirs

from spots.spots import coldspots, hotspots

from chains.chains import Chains
from mutation_strategies.mutation_strategies import MutationStrategies

from kmer_utilities.filtering_kmers_utilities import filter_by_coverage
from sample_reader.standard_samples import concatenate_kmer_freq_matrices
from mutability_diversity.calculate_mutability_diversity \
    import calculate_mutability_fr, \
           calculate_mutability_cdr, \
           calculate_mutability_full


def plot_mutability_boxplots(matrices, output_dir):
    hspots = set(hotspots())
    cspots = set(coldspots())

    def plot_by_strategy_chain(matrices, chain, strategy, output_dir):
        matrices = filter_by_coverage(
            matrices, 
            threshold_function=lambda x, axis: np.nanmin(x[:, [1, 3]], axis=axis),
            mean_function=np.nanmedian,
            coverage_threshold=1500)
        mut_fr = calculate_mutability_fr(matrices)
        mut_cdr = calculate_mutability_cdr(matrices)
        mut_full = calculate_mutability_full(matrices)

        def boxplots(mut, kmer_names, fig_name):
            med = np.nanmedian(mut, axis=1)
            order = np.argsort(med)[::-1]
            kmer_names = np.array(kmer_names)
            kmer_names = kmer_names[order]
            mut = mut[order, :]
            mut = pd.DataFrame(mut.T, columns=kmer_names)

            colors = np.array(['black'] * len(kmer_names))
            for i, kmer in enumerate(kmer_names):
                if kmer in hspots:
                    colors[i] = 'red'
                elif kmer in cspots:
                    colors[i] = 'blue'

            if not mut.empty:
                g = sns.boxplot(data=mut)
                g.set_xticklabels(g.get_xticklabels(), rotation=90)
                [t.set_color(i) for (i,t) in zip(colors, g.xaxis.get_ticklabels())]
                [t.set_facecolor(i) for (i,t) in zip(colors, g.artists)]
                g.set_ylim(0, 1)
                fig = g.get_figure()
                fig.set_size_inches(10, 7)
                fig.savefig(fig_name, format='pdf', dpi=150)
                plt.close()

        kmer_names = matrices.kmer_names
        fig_name_fr = "output_%s_%s_fr.pdf" % (strategy.name, chain.name)
        fig_name_cdr = "output_%s_%s_cdr.pdf" % (strategy.name, chain.name)
        fig_name_full = "output_%s_%s_full.pdf" % (strategy.name, chain.name)

        fig_name_fr = os.path.join(output_dir, fig_name_fr)
        fig_name_cdr = os.path.join(output_dir, fig_name_cdr)
        fig_name_full = os.path.join(output_dir, fig_name_full)

        boxplots(mut_fr, kmer_names, fig_name_fr)
        boxplots(mut_cdr, kmer_names, fig_name_cdr)
        boxplots(mut_full, kmer_names, fig_name_full)

    for strategy in MutationStrategies:
        for chain in Chains:
            print(strategy, chain)
            plot_by_strategy_chain(matrices[strategy][chain],
                                   chain, strategy, output_dir)


if __name__ == "__main__":
    input_config = read_config(parse_args().input)
    matrices = concatenate_kmer_freq_matrices(
        input_data=input_config.input_data,
        prefix_dir=input_config.prefix_dir,
        dir_data=input_config.kmer_model_estimating.kmer_matrices_dir,
        filename_fr=input_config.kmer_model_estimating.filename_fr,
        filename_cdr=input_config.kmer_model_estimating.filename_cdr)

    spots_boxplots = os.path.join(input_config.kmer_model_estimating.figures_dir,
                                  input_config.kmer_model_estimating.spots_boxplots)
    smart_makedirs(spots_boxplots)
    plot_mutability_boxplots(matrices, spots_boxplots)
