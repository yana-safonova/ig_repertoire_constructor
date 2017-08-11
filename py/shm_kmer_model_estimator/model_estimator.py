import os.path

import numpy as np
import pandas as pd

from special_utils.os_utils import smart_makedirs

from config.config import read_config
from config.parse_input_args import parse_args
from sample_reader.standard_samples import concatenate_kmer_freq_matrices
from shm_kmer_model_estimator.shm_kmer_model_estimator \
    import ShmKmerModelEstimator

from experiments.model_analysis import convergence_analysis
from experiments.spots_boxplots import plot_mutability_boxplots
from experiments.fr_cdr_comparison import compare_fr_cdr
from experiments.loci_comparison import compare_loci
from experiments.forward_backward_mutations import compare_forward_backward_mutations
from experiments.plot_beta_dir import plot_distr
from experiments.model_analysis import read_models
from experiments.ci_analysis import ci_len_analysis

from chains.chains import Chains
from shm_kmer_model.cab_shm_model import Region, CAB_SHM_Model
from mutation_strategies.mutation_strategies import MutationStrategies
from special_utils.model_utils import stddict


def output_models(models, output_directory):
    smart_makedirs(output_directory)
    for strategy in models:
        for chain in models[strategy]:
            model = models[strategy][chain]
            path = os.path.join(output_directory,
                                strategy.name + '_' + chain.name + '.csv')
            model.to_csv(path, na_rep='nan')


def wrapper_convergence_analysis(models, model_config):
    chains = [Chains.IGH, Chains.IGK, Chains.IGL]
    mut_str = [MutationStrategies.NoKNeighbours]
    outdir = os.path.join(model_config.outdir,
                          model_config.analysis_dir)
    outfile = os.path.join(outdir, model_config.model_convergence_analysis)
    smart_makedirs(outdir)
    convergence_analysis(models, verbose=False,
                         chains=chains,
                         strategies=mut_str).to_csv(outfile, sep=',')


def wrapper_plot_mutability_boxplots(matrices, model_config):
    spots_boxplots = os.path.join(model_config.figures_dir,
                                  model_config.spots_boxplots)
    smart_makedirs(spots_boxplots)
    plot_mutability_boxplots(matrices, spots_boxplots)


def make_figures_analysis_dirs(model_config, figures_exp_dir):
    figures_dir = os.path.join(model_config.figures_dir,
                               figures_exp_dir)
    smart_makedirs(figures_dir)

    outdir = os.path.join(model_config.outdir,
                          model_config.analysis_dir)
    smart_makedirs(outdir)
    return figures_dir, outdir


def wrapper_compare_fr_cdr(matrices, model_config):
    figures_dir, outdir = make_figures_analysis_dirs(
        model_config,
        model_config.fr_cdr_comparison_dir)
    outfile = os.path.join(outdir, model_config.fr_cdr_comparison_filename)
    compare_fr_cdr(matrices, figures_dir).to_csv(outfile, sep=',')


def wrapper_compare_loci(matrices, model_config):
    figures_dir, outdir = make_figures_analysis_dirs(
        model_config,
        model_config.loci_comparison_dir)
    outfile = os.path.join(outdir, model_config.loci_comparison_filename)
    compare_loci(matrices, figures_dir).to_csv(outfile, sep=',')


def wrapper_compare_forward_backward_mutations(matrices, model_config):
    figures_dir, outdir = make_figures_analysis_dirs(
        model_config,
        model_config.forward_backward_mutations_comparison_dir)
    outfile = os.path.join(
        outdir,
        model_config.forward_backward_mutations_comparison_filename)
    compare_forward_backward_mutations(matrices,
                                       figures_dir).to_csv(outfile, sep=',')


def wrapper_plot_distr(models, model_config):
    figures_dir = os.path.join(model_config.figures_dir,
                               model_config.model_distribution_figs)
    plot_distr(models, figures_dir)


def bootstrap(models, matrices, model_config):
    results, rand_matrices = stddict([]), stddict([])

    for b in xrange(model_config.bootstrap_iterations):
        print("%d / %d" % (b + 1, model_config.bootstrap_iterations))
        for strategy in matrices:
            for chain in matrices[strategy]:
                matrix = matrices[strategy][chain]
                rand_matrices[strategy][chain] = matrix.get_bootstrap_matrices()

        models = ShmKmerModelEstimator(model_config).estimate_models(rand_matrices)

        for strategy in matrices:
            for chain in matrices[strategy]:
                model = CAB_SHM_Model(dataset=models[strategy][chain],
                                      strategy=strategy, chain=chain)
                results[strategy][chain].append(model)
                # print(model.dataset.head())

        # for strategy in matrices:
        #     for chain in matrices[strategy]:
        #         print(results[strategy][chain][0].dataset.head())


    for strategy in matrices:
        for chain in matrices[strategy]:
            means_beta_fr, means_beta_cdr, means_beta_full = [], [], []
            means_dir1, means_dir2, means_dir3 = [], [], []
            for i in xrange(model_config.bootstrap_iterations):
                model = results[strategy][chain][i]
                means_beta_fr.append(model.all_expect_mut_probs(region=Region.FR))
                means_beta_cdr.append(model.all_expect_mut_probs(region=Region.CDR))
                means_beta_full.append(model.all_expect_mut_probs())
                all_dir_means = np.array(model.all_dir_means())
                means_dir1.append(all_dir_means[:,0])
                means_dir2.append(all_dir_means[:,1])
                means_dir3.append(all_dir_means[:,2])

            means_beta_fr = np.stack(means_beta_fr, axis=1)
            means_beta_cdr = np.stack(means_beta_cdr, axis=1)
            means_beta_full = np.stack(means_beta_full, axis=1)
            means_dir1 = np.stack(means_dir1, axis=1)
            means_dir2 = np.stack(means_dir2, axis=1)
            means_dir3 = np.stack(means_dir3, axis=1)

            percentile_beta_fr = np.nanpercentile(means_beta_fr, [5, 95], axis=1)
            percentile_beta_cdr = np.nanpercentile(means_beta_cdr, [5, 95], axis=1)
            percentile_beta_full = np.nanpercentile(means_beta_full, [5, 95], axis=1)
            percentile_dir1 = np.nanpercentile(means_dir1, [5, 95], axis=1)
            percentile_dir2 = np.nanpercentile(means_dir2, [5, 95], axis=1)
            percentile_dir3 = np.nanpercentile(means_dir3, [5, 95], axis=1)

            model = models[strategy][chain]

            def add_ci(df, what, name):
                df[name + "_left"] = pd.Series(what[0, :], index=df.index)
                df[name + "_right"] = pd.Series(what[1, :], index=df.index)
                return df

            model = add_ci(model, percentile_beta_fr, "percentile_beta_fr")
            model = add_ci(model, percentile_beta_cdr, "percentile_beta_cdr")
            model = add_ci(model, percentile_beta_full, "percentile_beta_full")
            model = add_ci(model, percentile_dir1, "percentile_dir1")
            model = add_ci(model, percentile_dir2, "percentile_dir2")
            model = add_ci(model, percentile_dir3, "percentile_dir3")
    return models


def main():
    parsed_args = parse_args()
    print("Reading input config from %s" % parsed_args.input)
    input_config = read_config(parsed_args.input)
    model_config = input_config.kmer_model_estimating
    print("Reading kmer matrices started")
    matrices = concatenate_kmer_freq_matrices(
        input_data=input_config.input_data,
        prefix_dir=input_config.prefix_dir,
        dir_data=model_config.kmer_matrices_dir,
        functionality=model_config.functionality,
        filename_fr=model_config.filename_fr,
        filename_cdr=model_config.filename_cdr)
    print("Reading kmer matrices ended")

    wrapper_compare_forward_backward_mutations(matrices, model_config)
    print("Estimating models started")
    models = ShmKmerModelEstimator(model_config).estimate_models(matrices)
    print("Estimating models finished")

    print("Bootstrap CIs for model params")
    models = bootstrap(models, matrices, model_config)

    print("Outputing models to %s" % model_config.outdir)
    output_models(models, model_config.outdir)

    print("Analysis of models' convergence")
    wrapper_convergence_analysis(models, model_config)

    print("Analysis of ci_len")
    ci_len_analysis(model_config)

    if not parsed_args.skip_analysis:
        print("Plotting mutability boxplots started")
        wrapper_plot_mutability_boxplots(matrices, model_config)
        print("Plotting mutability boxplots finished")

        print("FR and CDR comparison started")
        wrapper_compare_fr_cdr(matrices, model_config)
        print("FR and CDR comparison finished")

        print("Various loci comparison started")
        wrapper_compare_loci(matrices, model_config)
        print("Various loci comparison ended")

        print("Comparison of forward and backward mutations started")
        wrapper_compare_forward_backward_mutations(matrices, model_config)
        print("Comparison of forward and backward mutations finished")

        models = read_models(model_config.outdir)
        print("Model distribution plots are drawing...")
        wrapper_plot_distr(models, model_config)
        print("Finished drawing model distribution plots")


if __name__ == "__main__":
    main()
