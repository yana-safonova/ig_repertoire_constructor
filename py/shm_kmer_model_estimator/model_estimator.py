import os.path
import json

from special_utils.os_utils import smart_makedirs

from config.config import config, read_config
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

from chains.chains import Chains
from mutation_strategies.mutation_strategies import MutationStrategies


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
    models = ShmKmerModelEstimator().estimate_models(matrices)
    print("Estimating models finished")


    print("Outputing models to %s" % model_config.outdir)
    output_models(models, model_config.outdir)

    print("Analysis of models' convergence")
    wrapper_convergence_analysis(models, model_config)

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
