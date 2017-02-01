import os.path

from config.config import config
from sample_reader.standard_samples import concatenate_kmer_matrices_all_data
from shm_kmer_model_estimator.shm_kmer_model_estimator import ShmKmerModelEstimator


def output_models(models):
    try: 
        os.makedirs(config.output_directory)
    except OSError:
        if not os.path.isdir(config.output_directory):
            raise
    for strategy in models:
        for chain in models[strategy]:
            model = models[strategy][chain]
            path = os.path.join(config.output_directory,
                                strategy.name + '_' + chain.name + '.csv')
            model.dataset.to_csv(path, na_rep='nan')


def main():
    matrices = concatenate_kmer_matrices_all_data()
    models = ShmKmerModelEstimator().estimate_models(matrices)
    output_models(models)


if __name__ == "__main__":
    main()
