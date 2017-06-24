from __future__ import division

import os.path
from collections import defaultdict

import numpy as np
import pandas as pd
from statsmodels.stats.proportion import proportion_confint

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['figure.figsize'] = 15, 10

from special_utils.os_utils import smart_makedirs
from config.config import config, read_config
from config.parse_input_args import parse_args

from spots.spots import coldspots, hotspots
from sample_reader.standard_samples import concatenate_kmer_freq_matrices
from shm_kmer_model.cab_shm_model import Region
from shm_kmer_model.cab_shm_nonproductive_model import CAB_SHM_Nonproductive_Model


def sum_matrices(config, matrices):
    summed_matrices_dir = os.path.join(config.outdir, config.summed_matrices_dir)
    smart_makedirs(summed_matrices_dir)
    for k1 in matrices:
        for k2 in matrices[k1]:
            matrices[k1][k2].matrices = np.sum(matrices[k1][k2].matrices, axis=2)
            filename = os.path.join(summed_matrices_dir, "%s_%s.csv" % (k1.name, k2.name))
            np.savetxt(filename, matrices[k1][k2].matrices, delimiter=";", fmt="%d")
    return matrices


def estimate_params(config, matrices):
    def estimate_params__(config, matrix):
        fr_mut   =  matrix[:, 0] / (matrix[:, 0] + matrix[:, 1])
        cdr_mut  =  matrix[:, 2] / (matrix[:, 2] + matrix[:, 3])
        full_mut = (matrix[:, 0] + matrix[:, 2]) / np.sum(matrix[:, 0:4], axis=1)
        subst    = (matrix[:, 4:].T / np.sum(matrix[:, 4:], axis=1)).T

        def get_mut_ci(counts, nobs, alpha=config.ci_mut_alpha, method=config.ci_mut_method):
            res = proportion_confint(counts, nobs, alpha, method=method)
            return np.array([(res[0][i], res[1][i]) for i in xrange(len(res[0]))])

        fr_mut_ci   = get_mut_ci(matrix[:, 0], matrix[:, 0] + matrix[:, 1])
        cdr_mut_ci  = get_mut_ci(matrix[:, 2], matrix[:, 2] + matrix[:, 3])
        full_mut_ci = get_mut_ci(matrix[:, 0] + matrix[:, 2],
                                 np.sum(matrix[:, 0:4], axis=1))
        return fr_mut, fr_mut_ci, cdr_mut, cdr_mut_ci, full_mut, full_mut_ci, subst

    models = defaultdict(lambda: defaultdict(int))
    for k1 in matrices:
        for k2 in matrices[k1]:
            fr_mut, fr_mut_ci, cdr_mut, cdr_mut_ci, full_mut, full_mut_ci, subst = \
                estimate_params__(config, matrices[k1][k2].matrices)
            model_dataset = np.concatenate((fr_mut[:, np.newaxis],
                                            fr_mut_ci,
                                            cdr_mut[:, np.newaxis],
                                            cdr_mut_ci,
                                            full_mut[:, np.newaxis],
                                            full_mut_ci,
                                            subst), axis=1)
            models[k1][k2] = CAB_SHM_Nonproductive_Model(strategy=k1,
                                                         chain=k2,
                                                         dataset=model_dataset,
                                                         export=True)
    return models


def draw_muts(models, config):
    def draw_mut_strategy_chain(model, config, strategy, chain, mut_plots_dir):
        def draw_mut_region(model, config, muts, cis, kmer_ind, kmers, figname):
            muts = muts[kmer_ind]
            cis = cis[kmer_ind]

            order = np.argsort(muts)[::-1]
            muts = muts[order]
            kmer_ind = kmer_ind[order]
            kmers = kmers[order]
            cis = cis[order]
            muts = np.concatenate((muts[:, np.newaxis], cis), axis=1)

            muts = pd.DataFrame(muts.T, columns=kmers)
            print(figname)

            colors = np.array(['black'] * len(kmers))
            hspots = set(hotspots())
            cspots = set(coldspots())
            for i, kmer in enumerate(kmers):
                if kmer in hspots:
                    colors[i] = 'red'
                elif kmer in cspots:
                    colors[i] = 'blue'

            g = sns.pointplot(data=muts, markers='_', linestyles='None', scale=0.1, join=False)
            g.set_xticklabels(g.get_xticklabels(), rotation=90)
            [t.set_color(i) for (i,t) in zip(colors, g.xaxis.get_ticklabels())]
            [t.set_facecolor(i) for (i,t) in zip(colors, g.artists)]
            g.set_ylim(0, 1)
            fig = g.get_figure()
            fig.set_size_inches(10, 7)
            fig.savefig(figname, format='pdf', dpi=150)
            plt.close()

        def call_draw(region):
            kmer_ind, kmers = model.get_well_covered_kmers(region=region)
            draw_mut_region(model, config,
                            muts = model.all_expect_mut_prob(region=region),
                            cis = model.all_expect_mut_ci(region=region),
                            kmer_ind=kmer_ind,
                            kmers=kmers,
                            figname=os.path.join(mut_plots_dir, '%s_%s_%s.pdf' % (strategy.name, chain.name, region.name)))

        call_draw(Region.FR)
        call_draw(Region.CDR)
        call_draw(Region.ANY)



    mut_plots_dir = os.path.join(config.outdir, config.mut_plots)
    smart_makedirs(mut_plots_dir)
    for k1 in models:
        for k2 in models[k1]:
            draw_mut_strategy_chain(models[k1][k2], config, k1, k2, mut_plots_dir)


def main():
    parsed_args = parse_args()
    print("Reading input config from %s" % parsed_args.input)
    input_config = read_config(parsed_args.input)
    model_config = input_config.kmer_model_estimating
    nonproductive_model_config = input_config.kmer_nonproductive_model_estimating
    print("Reading kmer matrices started")
    matrices = concatenate_kmer_freq_matrices(
        input_data=input_config.input_data,
        prefix_dir=input_config.prefix_dir,
        dir_data=model_config.kmer_matrices_dir,
        functionality=nonproductive_model_config.functionality,
        filename_fr=model_config.filename_fr,
        filename_cdr=model_config.filename_cdr)
    print("Reading kmer matrices ended")

    print("Sum kmer matrices")
    matrices = sum_matrices(config=nonproductive_model_config, matrices=matrices)

    print("Estimating params")
    models = estimate_params(config=nonproductive_model_config,
                             matrices=matrices)

    print("Plotting Mutabilities with CIs")
    draw_muts(config=nonproductive_model_config, models=models)

if __name__ == "__main__":
    main()
