from enum import Enum

import numpy as np
import os
import pandas as pd

from abc import abstractmethod
import kmer_utilities.kmer_utilities as kmer_utilities
import mutation_strategies.mutation_strategies as mutation_strategies
import chains.chains as chains

from config.config import config


class Region(Enum):
    FR, CDR, ANY = range(3)


class CAB_SHM_Model(object):
    model_path = os.path.join('/', 'Sid', 'abzikadze', 'shm_models')

    def __init__(self, strategy, chain, kmer_len=5):
        self.kmer_len = kmer_len
        self.kmer_names = kmer_utilities.kmer_names(kmer_len)
        self.column_names = config.output_csv_header
        self.strategy = strategy 
        self.chain = chain
        path = os.path.join(self.model_path, '%s_%s.csv' % (strategy.name, chain.name))
        self.dataset = pd.read_csv(path, index_col=0)

    def __get_beta_params(self, region=Region.ANY):
        if region == Region.FR:
            beta_shape1, beta_shape2 = 'beta_FR_shape1', 'beta_FR_shape2'
        elif region == Region.CDR:
            beta_shape1, beta_shape2 = 'beta_CDR_shape1', 'beta_CDR_shape2'
        elif region == Region.ANY:
            beta_shape1, beta_shape2 = 'beta_FULL_shape1', 'beta_FULL_shape2'
        else:
             raise ValueError("Wrong region")
        return beta_shape1, beta_shape2

    def beta_params(self, kmer, region=Region.ANY):
        assert kmer in self.kmer_names
        p1, p2 = self.__get_beta_params(region)
        return self.dataset.loc[kmer][[p1, p2]]

    def all_beta_params(self, region=Region.ANY):
        p1, p2 = self.__get_beta_params(region)
        return self.dataset[[p1, p2]]

    def dir_params(self, kmer):
        assert kmer in self.kmer_names
        return self.dataset.loc[kmer][['dir_shape1',
                                       'dir_shape2',
                                       'dir_shape3']]

    def all_dir_params(self):
        return self.dataset[['dir_shape1',
                             'dir_shape2',
                             'dir_shape3']]

    def __get_start_point_params(self, region=Region.ANY):
        if region == Region.FR:
             p1 = 'start_point_beta_FR_shape1'
             p2 = 'start_point_beta_FR_shape2'
        elif region == Region.CDR:
             p1 = 'start_point_beta_CDR_shape1'
             p2 = 'start_point_beta_CDR_shape2'
        elif region == Region.ANY:
             p1 = 'start_point_beta_FULL_shape1'
             p2 = 'start_point_beta_FULL_shape2'
        else:
             raise ValueError("Wrong region")
        return p1, p2

    def start_point_beta_params(self, kmer, region=Region.ANY):
        assert kmer in self.kmer_names
        p1, p2 = self.__get_start_point_params(region)
        return self.dataset.loc[kmer][[p1, p2]]

    def all_start_point_beta_params(self, region):
        p1, p2 = self.__get_start_point_params(region)
        return self.dataset[[p1, p2]]

    def start_point_dir_params(self, kmer):
        assert kmer in self.kmer_names
        return self.dataset.loc[kmer][['start_point_dir_shape1',
                                       'start_point_dir_shape2',
                                       'start_point_dir_shape3']]

    def all_start_point_dir_params(self):
        return self.dataset[['start_point_dir_shape1',
                             'start_point_dir_shape2',
                             'start_point_dir_shape3']]

    def expect_mut_prob(self, kmer, region=Region.ANY):
        assert kmer in self.kmer_names
        p1, p2 = self.__get_beta_params(region)
        kmer_info = self.dataset.loc[kmer]
        result = kmer_info[p1] / \
            (kmer_info[p1] + kmer_info[p2])
        if not np.isnan(result):
            return result
        return kmer_info['start_point_' + p1]

    def all_expect_mut_probs(self, region=Region.ANY):
        result = []
        for kmer in self.kmer_names:
            result.append(self.expect_mut_prob(kmer, region))
        return np.array(result)

    def expect_subst_prob(self, kmer, subst_base):
        assert kmer in self.kmer_names

        bases = kmer_utilities.nucl_bases()
        assert subst_base in bases

        central_nucl = kmer[self.kmer_len // 2]
        assert central_nucl != subst_base

        bases.remove(central_nucl)
        mutation_index = bases.index(subst_base) + 1

        kmer_info = self.dataset.loc[kmer]
        result = kmer_info['dir_shape' + str(mutation_index)] / \
            (kmer_info['dir_shape1'] +
             kmer_info['dir_shape2'] +
             kmer_info['dir_shape3'])
        if not np.isnan(result):
            return result
        return kmer_info['start_point_dir_shape' + str(mutation_index)]
