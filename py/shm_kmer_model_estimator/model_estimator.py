import os.path
import json

from special_utils.os_utils import smart_makedirs

from config.config import config, read_config
from config.parse_input_args import parse_args
from sample_reader.standard_samples import concatenate_kmer_freq_matrices
from shm_kmer_model_estimator.shm_kmer_model_estimator \
    import ShmKmerModelEstimator

from experiments.model_analysis import coverage_kmers_utils_all_models
from experiments.spots_boxplots import plot_mutability_boxplots

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


def convergence_analysis(models, input_config):
    chains = [Chains.IGH, Chains.IGK, Chains.IGL]
    mut_str = [MutationStrategies.NoKNeighbours]
    model_config = input_config.kmer_model_estimating
    outdir = os.path.join(model_config.outdir,
                          model_config.analysis_dir)
    outfile = os.path.join(outdir, model_config.kmer_coverage_analysis)
    smart_makedirs(outdir)
    coverage_kmers_utils_all_models(models, verbose=False,
                                    chains=chains,
                                    strategies=mut_str).to_csv(outfile, sep=',')


def mutability_boxplots(matrices, input_config):
    spots_boxplots = os.path.join(input_config.kmer_model_estimating.figures_dir,
                                  input_config.kmer_model_estimating.spots_boxplots)
    smart_makedirs(spots_boxplots)
    plot_mutability_boxplots(matrices, spots_boxplots)


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
    mutability_boxplots(matrices, input_config)
    print("Plotting mutability boxplots finished")

    print("Estimating models started")
    models = ShmKmerModelEstimator().estimate_models(matrices)
    print("Estimating models finished")
    print("Outputing models to %s" % input_config.kmer_model_estimating.outdir)
    output_models(models, input_config.kmer_model_estimating.outdir)

    print("Analysis of models' convergence")
    convergence_analysis(models, input_config)


if __name__ == "__main__":
    main()
