#!/usr/bin/env python2

from __future__ import print_function

import numpy as np
import pandas as pd

from abc import abstractmethod
import kmer_utilities.kmer_utilities as kmer_utilities


class AbstractSHM_Model(object):
    def __init__(self, dataset, kmer_len=5):
        self.dataset = dataset
        self.kmer_len = kmer_len
        self.kmer_names = kmer_utilities.kmer_names(kmer_len)
        assert (self.kmer_names == self.dataset.index).all()

    def get_raw_dataset(self):
        return self.dataset

    def __str__(self):
        return self.dataset.__str__()

    def __repr__(self):
        return self.dataset.__repr__()

    @abstractmethod
    def expect_mut_prob(self, kmer):
        pass

    @abstractmethod
    def all_expect_mut_probs(self, method='all'):
        pass

    @abstractmethod
    def expect_subst_prob(self, kmer, subst_base):
        pass

    @abstractmethod
    def all_expect_subst_probs(self, method='all'):
        pass


class CAB_SHM_Model(AbstractSHM_Model):
    def __init__(self, dataset, kmer_len=5):
        super(CAB_SHM_Model, self).__init__(dataset, kmer_len)

    def beta_params(self, kmer):
        assert kmer in self.kmer_names
        return self.dataset.loc[kmer][['beta_shape1',
                                       'beta_shape2']]

    def all_beta_params(self):
        return self.dataset[['beta_shape1',
                             'beta_shape2']]

    def dir_params(self, kmer):
        assert kmer in self.kmer_names
        return self.dataset.loc[kmer][['dir_shape1',
                                       'dir_shape2',
                                       'dir_shape3']]

    def all_dir_params(self):
        return self.dataset[['dir_shape1',
                             'dir_shape2',
                             'dir_shape3']]

    def start_point_beta_params(self, kmer):
        assert kmer in self.kmer_names
        return self.dataset.loc[kmer][['start_point_beta_shape1',
                                       'start_point_beta_shape2']]

    def all_start_point_beta_params(self):
        return self.dataset[['start_point_beta_shape1',
                             'start_point_beta_shape2']]

    def start_point_dir_params(self, kmer):
        assert kmer in self.kmer_names
        return self.dataset.loc[kmer][['start_point_dir_shape1',
                                       'start_point_dir_shape2',
                                       'start_point_dir_shape3']]

    def all_start_point_dir_params(self):
        return self.dataset[['dir_shape1',
                             'dir_shape2',
                             'dir_shape3']]

    def expect_mut_prob(self, kmer):
        assert kmer in self.kmer_names
        kmer_info = self.dataset.loc[kmer]
        result = kmer_info['beta_shape1'] / \
                (kmer_info['beta_shape1'] + kmer_info['beta_shape2'])
        if not np.isnan(result):
            return result
        return kmer_info['start_point_beta_shape1']

    def all_expect_mut_probs(self):
        result = self.dataset['beta_shape1'] / \
                (self.dataset['beta_shape1'] + self.dataset['beta_shape2'])
        # TODO try to change NaN to start points
        return result

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
                (kmer_info['dir_shape1'] + \
                 kmer_info['dir_shape2'] + \
                 kmer_info['dir_shape3'])
        if not np.isnan(result):
            return result
        return kmer_info['start_point_dir_shape' + str(mutation_index)]


class YaleSHM_Model(AbstractSHM_Model):
    import os
    mutability_path = os.path.join('yale_model', 'Mutability.csv')
    substitution_path = os.path.join('yale_model', 'Substitution.csv')

    def __init__(self):
        self.mutability_dataset = pd.read_csv(self.mutability_path, sep=' ')
        self.mutability_dataset.set_index(self.mutability_dataset.Fivemer, inplace=True)
        self.mutability_dataset = self.mutability_dataset.drop('Fivemer', 1)

        self.substitution_dataset = pd.read_csv(self.substitution_path, sep=' ')
        self.substitution_dataset.set_index(self.substitution_dataset.Fivemer, inplace=True)
        self.substitution_dataset = self.substitution_dataset.drop('Fivemer', 1)

    def __str__(self):
        return self.mutability_dataset.__str__() + self.substitution_dataset.__str__()

    def __repr__(self):
        return self.mutability_dataset.__repr__() + self.substitution_dataset.__repr__()

    def expect_mut_prob(self, kmer):
        assert kmer in self.kmer_names
        kmer_info = self.mutability_dataset.loc[kmer]
        return kmer_info['Mutability']

    def __get_indexes_by_method(self, method):
        assert method == 'Measured' or method == 'Inferred'
        return self.mutability_dataset.Source == method

    def all_expect_mut_probs(self, method='all'):
        if method != 'all':
            indexes = self.__get_indexes_by_method(method)
            return self.mutability_dataset.Mutability.loc[indexes]
        return self.mutability_dataset

    def expect_subst_prob(self, kmer, subst_base):
        assert kmer in self.kmer_names

        bases = kmer_utilities.nucl_bases()
        assert subst_base in bases
        return self.substitution_dataset.loc[kmer][subst_base]

    def all_expect_subst_probs(self, method='all'):
        bases = kmer_utilities.nucl_bases()
        if method != 'all':
            indexes = self.__get_indexes_by_method(method)
            return self.substitution_dataset.loc[indexes][bases]
        return self.substitution_dataset[bases]
