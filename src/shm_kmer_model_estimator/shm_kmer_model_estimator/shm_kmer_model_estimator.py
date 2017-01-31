import pandas as pd
import numpy as np
from shm_kmer_likelihood_optimize.shm_kmer_likelihood_optimize \
    import ShmKmerLikelihoodOptimizator
from shm_kmer_likelihood.shm_kmer_likelihood import ShmKmerLikelihood
import kmer_utilities.kmer_utilities as kmer_utils
from shm_kmer_model.cab_shm_model import CAB_SHM_Model


class ShmKmerModelEstimator(object):
    def __init__(self, kmer_len=5):
        self.bases = kmer_utils.nucl_bases()
        self.n_nucl = len(self.bases)
        self.kmer_len = kmer_len
        self.half_kmer_len = kmer_len // 2
        self.kmer_names = kmer_utils.kmer_names(kmer_len)
        self.column_names = ['beta_FR_shape1', 'beta_FR_shape2',
                             'beta_CDR_shape1', 'beta_CDR_shape2',
                             'dir_shape1', 'dir_shape2', 'dir_shape3',
                             'success_optim_beta_FR',
                             'success_optim_beta_CDR',
                             'success_optim_dir',
                             'start_point_beta_FR_shape1',
                             'start_point_beta_FR_shape2',
                             'start_point_beta_CDR_shape1',
                             'start_point_beta_CDR_shape2',
                             'start_point_dir_shape1',
                             'start_point_dir_shape2',
                             'start_point_dir_shape3']
        self.central_nucl_indexes = kmer_utils.central_nucl_indexes(kmer_len)
        self.n_param = len(self.column_names)

    def __estimate_backend(self, kmer_matrices):
        kmers = kmer_utils.kmer_names()
        results = np.empty((len(kmers), self.n_param))

        for i, kmer in enumerate(kmers):
            if i % 250 == 0:
                print(i)
            lkho = ShmKmerLikelihood(kmer_matrices[kmer])
            res = ShmKmerLikelihoodOptimizator(lkho).maximize()
            fr = res['optim_res_beta_fr']
            cdr = res['optim_res_beta_cdr']
            direchlet = res['optim_res_dir']
            success = np.array([fr.success, cdr.success, direchlet.success])
            start_fr = res['start_beta_fr']
            start_cdr = res['start_beta_cdr']
            start_dir = res['start_dir']
            starts = np.concatenate([start_fr, start_cdr, start_dir])
            results[i, ] = \
                np.concatenate((fr.x, cdr.x, direchlet.x,
                                success,
                                starts))
        results = pd.DataFrame(data=results,
                               index=self.kmer_names,
                               columns=self.column_names)
        return CAB_SHM_Model(results)

    def estimate_models(self, dict_kmer_matrices,
                        strategies=None, chains=None):
        if strategies is None:
            strategies = dict_kmer_matrices.keys()
        if chains is None:
            chains = dict_kmer_matrices[strategies[0]].keys()
        models = dict.fromkeys(strategies)
        for strategy in strategies:
            models[strategy] = dict.fromkeys(chains)

        for strategy in strategies:
            for chain in chains:
                print("%s, %s" % (strategy.name, chain.name))
                kmer_matrices = dict_kmer_matrices[strategy][chain]
                models[strategy][chain] = self.__estimate_backend(kmer_matrices)
        return models
