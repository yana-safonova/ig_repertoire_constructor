import pandas as pd
import numpy as np

from joblib import Parallel, delayed
from joblib.pool import has_shareable_memory

from config.config import config

from shm_kmer_likelihood_optimize.shm_kmer_likelihood_optimize \
    import ShmKmerLikelihoodOptimizator
from shm_kmer_likelihood.shm_kmer_likelihood import ShmKmerLikelihood
import kmer_utilities.kmer_utilities as kmer_utils
from chains.chains import Chains

np.seterr(divide='ignore', invalid='ignore')


class ShmKmerModelEstimator(object):
    def __init__(self, model_config, kmer_len=5):
        self.bases = kmer_utils.nucl_bases()
        self.n_nucl = len(self.bases)
        self.kmer_len = kmer_len
        self.half_kmer_len = kmer_len // 2
        self.kmer_names = kmer_utils.kmer_names(kmer_len)
        self.column_names = config.output_csv_header
        self.central_nucl_indexes = kmer_utils.central_nucl_indexes(kmer_len)
        self.n_param = len(self.column_names)
        for chain in Chains:
            if chain == Chains.IG:
                continue
        self.min_mut_coverage = model_config.min_mut_coverage
        self.min_subst_coverage = model_config.min_subst_coverage


    def refine_model(self, results, matrices, chain):
        ind_need_mut_inference, ind_need_subst_inference = [], []
        for i, kmer in enumerate(kmer_utils.kmer_names()):
            if np.isnan(results.loc[kmer, "start_point_beta_FULL_shape1"]) or\
               np.isnan(results.loc[kmer, "start_point_beta_FULL_shape2"]):
                ind_need_mut_inference.append(i)

            if np.isnan(results.loc[kmer, "start_point_dir_shape1"]) or\
               np.isnan(results.loc[kmer, "start_point_dir_shape2"]) or\
               np.isnan(results.loc[kmer, "start_point_dir_shape3"]):
                ind_need_subst_inference.append(i)

        def update_res(i, res_ind, estimated_name):
            kmer = kmer_utils.kmer_names()[i]
            surr_kmers, surr_ind = \
                kmer_utils.surround_kmers(i=i, bad_ind=ind_need_mut_inference)
            surr_ind = np.array(surr_ind)
            means = np.array(results.iloc[surr_ind, res_ind])
            sums = np.nansum(means, axis=1)
            means = (means.T / sums).T
            mean = np.nanmean(means, axis=0)
            results.iloc[i, res_ind] = mean
            results.loc[kmer, estimated_name] = 0

        for i in ind_need_mut_inference:
            update_res(i, [4, 5], "beta_estimated")

        for i in ind_need_subst_inference:
            update_res(i, [6, 7, 8], "dir_estimated")

        return results

    def estimate_backend(self, kmer_matrices, strategy, chain):
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
                                starts,
                                [1], [1]))

        results = pd.DataFrame(data=results,
                               index=self.kmer_names,
                               columns=self.column_names)
        # results = self.refine_model(results, kmer_matrices, chain)
        # return CAB_SHM_Model(results)
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
            (models, strategy, chain, self, kmer_matrices)
            for strategy in strategies
            for chain in chains
        )
        for strategy, chain, model in r:
            models[strategy][chain] = model
        # for strategy in strategies:
        #     for chain in chains:
        #         models[strategy][chain] = \
        #             self.estimate_backend(kmer_matrices[strategy][chain],
        #                                   strategy, chain)
        return models


def joblib_est(models, strategy, chain, estimator, kmer_matrices):
    return strategy, chain, \
           estimator.estimate_backend(kmer_matrices[strategy][chain],
                                      strategy, chain)
