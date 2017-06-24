import os.path

import pandas as pd
import numpy as np
from special_utils.os_utils import smart_makedirs
import kmer_utilities.kmer_utilities as kmer_utilities
from cab_shm_model import Region

from config.config import config

class CAB_SHM_Nonproductive_Model(object):
    model_path = os.path.join('/', 'Sid', 'abzikadze', 'shm_models', 'nonproductive_model')

    def __init__(self, strategy, chain, kmer_len=5, dataset=None, export=False):
        self.kmer_len = kmer_len
        self.kmer_names = kmer_utilities.kmer_names(kmer_len)
        self.column_names = config.output_csv_header_nonproductive
        self.strategy = strategy 
        self.chain = chain
        self.path = os.path.join(self.model_path, '%s_%s.csv' % (strategy.name, chain.name))

        if dataset is None:
            self.dataset = pd.read_csv(path, index_col=0)
        else:
            self.dataset = pd.DataFrame(dataset,
                                        columns=self.column_names,
                                        index=self.kmer_names)
            if export:
                smart_makedirs(self.model_path)
                self.dataset.to_csv(self.path, sep=";")

        self.orig_matrix = np.genfromtxt(os.path.join(self.model_path, 'summed_matrices',
                                         '%s_%s.csv' % (strategy.name, chain.name)), delimiter=';')

    def get_mut_column__(self, region):
        if region == Region.FR:
            return "FR_mut"
        elif region == Region.CDR:
            return "CDR_mut"
        return "FULL_mut"

    def expect_mut_prob(self, kmer, region=Region.ANY):
        assert kmer in self.kmer_names
        return self.dataset.loc[kmer][self.get_mut_column__(region)]

    def all_expect_mut_prob(self, region=Region.ANY):
        return self.dataset[self.get_mut_column__(region)]

    def expect_subst_prob(self, kmer, subst_base):
        assert kmer in self.kmer_names

        bases = kmer_utilities.nucl_bases()
        assert subst_base in bases

        central_nucl = kmer[self.kmer_len // 2]
        assert central_nucl != subst_base

        bases.remove(central_nucl)
        mutation_index = bases.index(subst_base) + 1

        return self.dataset.loc[kmer]["subst_" + str(mutation_index)]

    def expect_mut_ci(self, kmer, region=Region.ANY):
        assert kmer in self.kmer_names
        col_name = self.get_mut_column__(region) 
        return (self.dataset.loc[kmer][col_name + "_ci_left"],
                self.dataset.loc[kmer][col_name + "_ci_right"])

    def all_expect_mut_ci(self, region=Region.ANY):
        col_name = self.get_mut_column__(region) 
        return np.concatenate((self.dataset[col_name + "_ci_left"][:, np.newaxis],
            self.dataset[col_name + "_ci_right"][:, np.newaxis]), axis=1)

    def get_well_covered_kmers(self, threshold=100, region=Region.ANY):
        if region == Region.FR:
            sel_cols = self.orig_matrix[:, [0, 1]]
        elif region == Region.CDR:
            sel_cols = self.orig_matrix[:, [2, 3]]
        else:
            sel_cols = np.concatenate((np.sum(self.orig_matrix[:, [0, 2]], axis=1)[:, np.newaxis],
                                       np.sum(self.orig_matrix[:, [1, 3]], axis=1)[:, np.newaxis]), axis=1)

        good_ind = np.min(sel_cols, axis=1) > threshold
        return good_ind, np.array(self.kmer_names)[good_ind]
