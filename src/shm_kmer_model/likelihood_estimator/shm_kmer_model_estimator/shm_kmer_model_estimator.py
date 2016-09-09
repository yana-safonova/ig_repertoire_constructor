import pandas as pd
import numpy as np
from shm_kmer_likelihood_optimize.shm_kmer_likelihood_optimize \
    import ShmKmerLikelihoodOptimizator
from shm_kmer_likelihood.shm_kmer_likelihood import ShmKmerLikelihood
import kmer_utilities.kmer_utilities as kmer_utils
import shm_kmer_model.shm_kmer_model as shm_kmer_model


class ShmKmerModelEstimator:
    def __init__(self, kmer_len=5):
        self.bases = kmer_utils.nucl_bases()
        self.n_nucl = len(self.bases)
        self.kmer_len = kmer_len
        self.half_kmer_len = kmer_len // 2
        self.kmer_names = kmer_utils.kmer_names(kmer_len)
        self.column_names = ['beta_shape1', 'beta_shape2',
                             'dir_shape1', 'dir_shape2', 'dir_shape3',
                             'success_optim_beta', 'success_optim_dir',
                             'start_point_beta_shape1',
                             'start_point_beta_shape2',
                             'start_point_dir_shape1',
                             'start_point_dir_shape2',
                             'start_point_dir_shape3']
        self.central_nucl_indexes = kmer_utils.central_nucl_indexes(kmer_len)
        self.n_param = len(self.column_names)

    def estimate_model(self, tuple_datasets):
        all_samples = np.concatenate(tuple_datasets, axis=0)
        n_datasets, n_kmer, n_nucl = all_samples.shape

        results = np.empty((n_kmer, self.n_param))
        for i in xrange(n_kmer):
            lkho = ShmKmerLikelihood(all_samples[:, i, :],
                                     self.central_nucl_indexes[i])
            start_p_beta, start_p_dir, optim_res_beta, optim_res_dir = \
                ShmKmerLikelihoodOptimizator(lkho).maximize()
            results[i, ] = np.concatenate((optim_res_beta.x, optim_res_dir.x,
                                          np.array([optim_res_beta.success,
                                                    optim_res_dir.success]),
                                          start_p_beta, start_p_dir))
        return shm_kmer_model.CAB_SHM_Model(pd.DataFrame(data=results,
                                            index=self.kmer_names,
                                            columns=self.column_names))

    def estimate_models_of_one_type(self, tuple_datasets,
                                    strategies=None, chains=None):
        if strategies is None:
            strategies = tuple_datasets[0].keys()
        if chains is None:
            chains = tuple_datasets[0][strategies[0]].keys()
        models = dict.fromkeys(strategies)
        for strategy in strategies:
            models[strategy] = dict.fromkeys(chains)

        for strategy in strategies:
            for chain in chains:
                print("%s: %s" % (strategy, chain))
                datasets = tuple(dataset[strategy][chain]
                                 for dataset in tuple_datasets)
                est_model = self.estimate_model(datasets)
                models[strategy][chain] = est_model
        return models
