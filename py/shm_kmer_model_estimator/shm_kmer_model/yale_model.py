from enum import Enum

import numpy as np
import pandas as pd

import kmer_utilities.kmer_utilities as kmer_utilities


class YaleMethod(Enum):
    all, Measured, Inferred = xrange(3)


class YaleSHM_Model(object):
    import os
    curr_path = os.path.dirname(os.path.abspath(__file__))
    mut_path = os.path.join(curr_path, '..', 'yale_model', 'Mutability.csv')
    subst_path = os.path.join(curr_path, '..', 'yale_model', 'Substitution.csv')

    def __init__(self):
        self.kmer_len = 5
        self.kmer_names = kmer_utilities.kmer_names(self.kmer_len)
        self.mut = pd.read_csv(self.mut_path, sep=' ')
        self.mut.set_index(self.mut.Fivemer, inplace=True)
        self.mut = self.mut.drop('Fivemer', 1)
        self.mut['Mutability'] /= self.mut['Mutability'].max()

        self.subst = pd.read_csv(self.subst_path, sep=' ')
        self.subst.set_index(self.subst.Fivemer, inplace=True)
        self.subst = self.subst.drop('Fivemer', 1)

    def __str__(self):
        return self.mut.__str__() + \
               self.subst.__str__()

    def __repr__(self):
        return self.mut.__repr__() + \
               self.subst.__repr__()

    def expect_mut_prob(self, kmer):
        assert kmer in self.kmer_names
        kmer_info = self.mut.loc[kmer]
        return kmer_info['Mutability']

    def __get_indexes_by_method(self, method):
        assert method == YaleMethod.Measured or method == YaleMethod.Inferred
        return self.mut.Source == method.name

    def all_expect_mut_probs(self, method=YaleMethod.all):
        if method != YaleMethod.all:
            indexes = self.__get_indexes_by_method(method)
            return self.mut.Mutability.loc[indexes]
        return self.mut

    def expect_subst_prob(self, kmer, subst_base):
        assert kmer in self.kmer_names

        bases = kmer_utilities.nucl_bases()
        assert subst_base in bases
        return self.subst.loc[kmer][subst_base]

    def all_expect_subst_probs(self, method=YaleMethod.all):
        bases = kmer_utilities.nucl_bases()
        if method != YaleMethod.all:
            indexes = self.__get_indexes_by_method(method)
            return self.subst.loc[indexes][bases]
        return self.subst[bases]
