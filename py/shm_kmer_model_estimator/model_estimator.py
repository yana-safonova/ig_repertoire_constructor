import os.path
import json

from config.config import config, read_config
from config.parse_input_args import parse_args
from sample_reader.standard_samples import concatenate_kmer_freq_matrices
from shm_kmer_model_estimator.shm_kmer_model_estimator \
    import ShmKmerModelEstimator


def output_models(models, output_directory):
    try: 
        os.makedirs(output_directory)
    except OSError:
        if not os.path.isdir(output_directory):
            raise
    for strategy in models:
        for chain in models[strategy]:
            model = models[strategy][chain]
            path = os.path.join(output_directory,
                                strategy.name + '_' + chain.name + '.csv')
            model.dataset.to_csv(path, na_rep='nan')


def main():
    input_config = read_config(parse_args().input)
    matrices = concatenate_kmer_freq_matrices(
        input_data=input_config.input_data,
        prefix_dir=input_config.prefix_dir,
        dir_data=input_config.kmer_model_estimating.kmer_matrices_dir,
        filename_fr=input_config.kmer_model_estimating.filename_fr,
        filename_cdr=input_config.kmer_model_estimating.filename_cdr)
    models = ShmKmerModelEstimator().estimate_models(matrices)
    output_models(models, input_config.kmer_model_estimating.outdir)


if __name__ == "__main__":
    main()
