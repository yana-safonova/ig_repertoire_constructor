""" This module contains standard samples from AbVitro and age datasets."""

import pandas as pd

from config.config import config
from kmer_matrices_reader import KmerMatricesReader
from kmer_matrices.kmer_matrices import KmerMatrices
from chains.chains import Chains
from mutation_strategies.mutation_strategies import MutationStrategies


def concatenate_kmer_freq_matrices(
        input_data,
        prefix_dir, dir_data,
        filename_fr, filename_cdr,
        functionality):
    kmer_matrices_reader = \
        KmerMatricesReader(prefix_dir=prefix_dir,
                           dir_data=dir_data,
                           filename_fr=filename_fr,
                           filename_cdr=filename_cdr)
    matrices = []
    for x in input_data:
        matrix = kmer_matrices_reader.read(
            root_dir=x.prefix_dir,
            functionality=functionality,
            ignore_indiv_number=x.bad_datasets)
        matrices.append(matrix)
    dict_matrices = dict.fromkeys(MutationStrategies)
    for strategy in MutationStrategies:
        matrices_str = \
            [matrix[strategy] for matrix in matrices]
        dict_matrices[strategy] = dict.fromkeys(Chains)
        for chain in Chains:
            matrices_chain = \
                [matrix[chain] for matrix in matrices_str if chain in matrix]
            dict_matrices[strategy][chain] = \
                KmerMatrices.FromKmerMatricesList(matrices_chain)
    return dict_matrices
