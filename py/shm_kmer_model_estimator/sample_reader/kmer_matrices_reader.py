import numpy as np
import pandas as pd

from config.config import config
from sample_reader import read_samples
from kmer_matrices.kmer_matrices import KmerMatrices


def read_matrix(x):
    pandas_array = pd.read_csv(x, sep=';', header=0, index_col=0)
    return np.array(pandas_array)


class KmerMatricesReader:
    def __init__(self, dir_data=config.kmer_matrices_dir,
                 filename_fr=config.filename_fr,
                 filename_cdr=config.filename_cdr,
                 read_func=read_matrix):
        self.dir_data = dir_data
        self.filename_fr = filename_fr
        self.filename_cdr = filename_cdr
        self.read_func = read_func 

    def __post_process(self, fr_matrices, cdr_matrices):
        matrices = dict.fromkeys(fr_matrices)
        for strategy in fr_matrices:
            matrices[strategy] = dict.fromkeys(fr_matrices[strategy])
            for chain_type in fr_matrices[strategy]:
                fr = fr_matrices[strategy][chain_type]
                cdr = cdr_matrices[strategy][chain_type]
                matrices[strategy][chain_type] = KmerMatrices(fr, cdr)
        return matrices

    def read(self, root_dir, ignore_indiv_number=[],
             prefix_dir=config.data_dir):
        def read_matrix(filename_data):
            return read_samples(filename_data=filename_data,
                                dir_data=self.dir_data,
                                read_func=self.read_func,
                                root_dir=root_dir,
                                ignore_indiv_number=ignore_indiv_number,
                                prefix_dir=prefix_dir)

        fr_matrices = read_matrix(self.filename_fr)
        cdr_matrices = read_matrix(self.filename_cdr)


        matrices = self.__post_process(fr_matrices, cdr_matrices)
        return matrices
