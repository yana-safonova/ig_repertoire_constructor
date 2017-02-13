import pandas as pd
import numpy as np

from joblib import Parallel, delayed
from joblib.pool import has_shareable_memory

from shm_kmer_likelihood_optimize.shm_kmer_likelihood_optimize \
    import ShmKmerLikelihoodOptimizator
from shm_kmer_likelihood.shm_kmer_likelihood import ShmKmerLikelihood
import kmer_utilities.kmer_utilities as kmer_utils
from shm_kmer_model.cab_shm_model import CAB_SHM_Model

from config.config import config

class ShmKmerModelEstimator(object):
    def __init__(self, kmer_len=5):
        self.bases = kmer_utils.nucl_bases()
        self.n_nucl = len(self.bases)
        self.kmer_len = kmer_len
        self.half_kmer_len = kmer_len // 2
        self.kmer_names = kmer_utils.kmer_names(kmer_len)
        self.column_names = config.output_csv_header
        self.central_nucl_indexes = kmer_utils.central_nucl_indexes(kmer_len)
        self.n_param = len(self.column_names)

    def estimate_backend(self, kmer_matrices):
        kmers = kmer_utils.kmer_names()
        results = np.empty((len(kmers), self.n_param))

        for i, kmer in enumerate(kmers):
            lkho = ShmKmerLikelihood(kmer_matrices[kmer])
            res = ShmKmerLikelihoodOptimizator(lkho).maximize()

            start_fr = res['start_beta_fr']
            start_cdr = res['start_beta_cdr']
            start_full = res['start_beta_full']
            start_dir = res['start_dir']
            starts = np.concatenate([start_fr, start_cdr,
                                     start_full, start_dir])

            fr = res['optim_res_beta_fr']
            cdr = res['optim_res_beta_cdr']
            full = res['optim_res_beta_full']
            direchlet = res['optim_res_dir']

            def x_success(start, max_obj):
                success = max_obj.success
                if np.any(np.isnan(max_obj.x)):
                    success = 0
                return success, max_obj.x

            fr_success, fr_x = x_success(start_fr, fr)
            cdr_success, cdr_x = x_success(start_cdr, cdr)
            full_success, full_x = x_success(start_full, full)
            dir_success, dir_x = x_success(start_dir, direchlet)

            success = np.array([fr_success, cdr_success,
                                full_success, dir_success])

            results[i, ] = \
                np.concatenate((fr_x, cdr_x, full_x, dir_x,
                                success,
                                starts))
        results = pd.DataFrame(data=results,
                               index=self.kmer_names,
                               columns=self.column_names)
        #return CAB_SHM_Model(results)
        return results

    def estimate_models(self, kmer_matrices,
                        strategies=None, chains=None):
        if strategies is None:
            strategies = kmer_matrices.keys()
        if chains is None:
            chains = kmer_matrices[strategies[0]].keys()
        models = dict.fromkeys(strategies)
        for strategy in strategies:
            models[strategy] = dict.fromkeys(chains)

        n_models = len(strategies) * len(chains)
        r = Parallel(n_jobs=n_models, max_nbytes=1e10)(
            delayed(joblib_est, has_shareable_memory)
            (models, strategy, chain, self, kmer_matrices) \
            for strategy in strategies \
            for chain in chains
        )
        for strategy, chain, model in r:
            models[strategy][chain] = model
        # for strategy in strategies:
        #     for chain in chains:
        #         models[strategy][chain] = \
        #             self.estimate_backend(kmer_matrices[strategy][chain])
        return models

def joblib_est(models, strategy, chain, estimator, kmer_matrices):
    return strategy, chain, \
           estimator.estimate_backend(kmer_matrices[strategy][chain])
