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


def wrapper_convergence_analysis(models, input_config):
    chains = [Chains.IGH, Chains.IGK, Chains.IGL]
    mut_str = [MutationStrategies.NoKNeighbours]
    model_config = input_config.kmer_model_estimating
    outdir = os.path.join(model_config.outdir,
                          model_config.analysis_dir)
    outfile = os.path.join(outdir, model_config.model_convergence_analysis)
    smart_makedirs(outdir)
    convergence_analysis(models, verbose=False,
                         chains=chains,
                         strategies=mut_str).to_csv(outfile, sep=',')


def wrapper_plot_mutability_boxplots(matrices, input_config):
    spots_boxplots = os.path.join(input_config.kmer_model_estimating.figures_dir,
                                  input_config.kmer_model_estimating.spots_boxplots)
    smart_makedirs(spots_boxplots)
    plot_mutability_boxplots(matrices, spots_boxplots)


def wrapper_compare_fr_cdr(matrices, input_config):
    model_config = input_config.kmer_model_estimating
    figures_dir = os.path.join(model_config.figures_dir,
                               model_config.fr_cdr_comparison_dir)
    smart_makedirs(figures_dir)

    outdir = os.path.join(model_config.outdir,
                          model_config.analysis_dir)
    outfile = os.path.join(outdir, model_config.fr_cdr_comparison_filename)
    smart_makedirs(outdir)
    compare_fr_cdr(matrices, figures_dir).to_csv(outfile, sep=',')


def wrapper_compare_loci(matrices, input_config):
    model_config = input_config.kmer_model_estimating
    figures_dir = os.path.join(model_config.figures_dir,
                               model_config.loci_comparison_dir)
    smart_makedirs(figures_dir)
    outdir = os.path.join(model_config.outdir,
                          model_config.analysis_dir)
    outfile = os.path.join(outdir, model_config.loci_comparison_filename)
    smart_makedirs(outdir)
    compare_loci(matrices, figures_dir).to_csv(outfile, sep=',')


def main():
    print("Reading input config from %s" % parse_args().input)
    input_config = read_config(parse_args().input)
    print("Reading kmer matrices started")
    matrices = concatenate_kmer_freq_matrices(
        input_data=input_config.input_data,
        prefix_dir=input_config.prefix_dir,
        dir_data=input_config.kmer_model_estimating.kmer_matrices_dir,
        filename_fr=input_config.kmer_model_estimating.filename_fr,
        filename_cdr=input_config.kmer_model_estimating.filename_cdr)
    print("Reading kmer matrices ended")
    print("Plotting mutability boxplots started")
    wrapper_plot_mutability_boxplots(matrices, input_config)
    print("Plotting mutability boxplots finished")

    print("FR and CDR comparison started")
    wrapper_compare_fr_cdr(matrices, input_config)
    print("FR and CDR comparison finished")

    print("Various loci comparison started")
    wrapper_compare_loci(matrices, input_config)
    print("Various loci comparison ended")

    print("Estimating models started")
    models = ShmKmerModelEstimator().estimate_models(matrices)
    print("Estimating models finished")
    print("Outputing models to %s" % input_config.kmer_model_estimating.outdir)
    output_models(models, input_config.kmer_model_estimating.outdir)

    print("Analysis of models' convergence")
    wrapper_convergence_analysis(models, input_config)


if __name__ == "__main__":
    main()
