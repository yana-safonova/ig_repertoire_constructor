import os.path
import argparse
import json

from config.config import config
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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        help="Config with input files",
                        required=True)
    return parser.parse_args()


def read_input_config():
    args = parse_args()
    with open(args.input, "r") as handle:
        config = json.load(handle,
                           object_hook=lambda d: argparse.Namespace(**d))
    return config.prefix_dir, config.input_data, config.kmer_model_estimating


def main():
    prefix_dir, input_data, model_est_params = read_input_config()
    matrices = concatenate_kmer_freq_matrices(
        input_data=input_data,
        prefix_dir=prefix_dir,
        dir_data=model_est_params.kmer_matrices_dir,
        filename_fr=model_est_params.filename_fr,
        filename_cdr=model_est_params.filename_cdr)
    models = ShmKmerModelEstimator().estimate_models(matrices)
    output_models(models, model_est_params.outdir)


if __name__ == "__main__":
    main()
