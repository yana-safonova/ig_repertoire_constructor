from enum import Enum

import numpy as np
import pandas as pd

from abc import abstractmethod
import kmer_utilities.kmer_utilities as kmer_utilities


class ModelMode(Enum):
    Mutation, Substitution, Both = range(3)

class Region(Enum):
    FR, CDR = range(2)


class CAB_SHM_Model(object):
    def __init__(self, dataset, kmer_len=5):
        self.kmer_len = kmer_len
        self.kmer_names = kmer_utilities.kmer_names(kmer_len)
        self.dataset = dataset.copy()

    def dump(self, filename):
        import pickle
        try:
            pickle.dump(self, filename)
        except:
            raise ValueError('Wrong filename')

    def __get_beta_params(self, region):
        if region == Region.FR:
            beta_shape1, beta_shape2 = 'beta_FR_shape1', 'beta_FR_shape2'
        elif region == Region.CDR:
            beta_shape1, beta_shape2 = 'beta_CDR_shape1', 'beta_CDR_shape2'
        return beta_shape1, beta_shape2

    def beta_params(self, region, kmer):
        assert kmer in self.kmer_names
        p1, p2 = self.__get_beta_params(region)
        return self.dataset.loc[kmer][[p1, p2]]

    def all_beta_params(self, region):
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

    def __get_start_point_params(self, region):
        if region == Region.FR:
             p1 = 'start_point_beta_FR_shape1'
             p2 = 'start_point_beta_FR_shape2'
        elif region == Region.CDR:
             p1 = 'start_point_beta_CDR_shape1'
             p2 = 'start_point_beta_CDR_shape2'
        else:
            raise ValueError("Wrong region")
        return p1, p2

    def start_point_beta_params(self, region, kmer):
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
        return self.dataset[['dir_shape1',
                             'dir_shape2',
                             'dir_shape3']]

    # def expect_mut_prob(self, kmer, region):
    #     assert kmer in self.kmer_names
    #     kmer_info = self.dataset.loc[kmer]
    #     result = kmer_info['beta_shape1'] / \
    #         (kmer_info['beta_shape1'] + kmer_info['beta_shape2'])
    #     if not np.isnan(result):
    #         return result
    #     return kmer_info['start_point_beta_shape1']

    # def all_expect_mut_probs(self):
    #     result = self.dataset['beta_shape1'] / \
    #             (self.dataset['beta_shape1'] + self.dataset['beta_shape2'])
    #     # TODO try to change NaN to start points
    #     return result

    # def expect_subst_prob(self, kmer, subst_base):
    #     assert kmer in self.kmer_names

    #     bases = kmer_utilities.nucl_bases()
    #     assert subst_base in bases

    #     central_nucl = kmer[self.kmer_len // 2]
    #     assert central_nucl != subst_base

    #     bases.remove(central_nucl)
    #     mutation_index = bases.index(subst_base) + 1

    #     kmer_info = self.dataset.loc[kmer]
    #     result = kmer_info['dir_shape' + str(mutation_index)] / \
    #         (kmer_info['dir_shape1'] +
    #          kmer_info['dir_shape2'] +
    #          kmer_info['dir_shape3'])
    #     if not np.isnan(result):
    #         return result
    #     return kmer_info['start_point_dir_shape' + str(mutation_index)]
