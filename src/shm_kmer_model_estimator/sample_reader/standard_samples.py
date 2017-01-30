""" This module contains standard samples from AbVitro and age datasets."""

import pandas as pd

from kmer_matrices_reader import KmerMatricesReader
from kmer_matrices.kmer_matrices import KmerMatrices
from chains.chains import Chains
from mutation_strategies.mutation_strategies import MutationStrategies


def concatenate_kmer_freq_matrices(
        input_data,
        kmer_matrices_reader=KmerMatricesReader()):
    matrices = [kmer_matrices_reader.read(*x) for x in input_data]
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


def concatenate_kmer_matrices_all_data(
        kmer_matrices_reader=KmerMatricesReader()):
    input_data = [("AbVitro/flu_time_course/FV/", ['25']),
                  ("AbVitro/flu_time_course/GMC/", ['8']),
                  ("AbVitro/flu_time_course/IDO/", []),
                  ("AbVitro/paired", []),
                  ("age", [])]
    result = concatenate_kmer_freq_matrices(
        input_data=input_data,
        kmer_matrices_reader=kmer_matrices_reader)
    return result
