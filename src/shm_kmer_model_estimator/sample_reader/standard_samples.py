""" This module contains standard samples from AbVitro and age datasets.
Func concatenate_samples creates a Panel that is as full as possible. """

from kmer_freq_matrix_reader import KmerFreqMatrixReader

import pandas as pd


def concatenate_kmer_freq_matrices(
        input_data,
        kmer_freq_matrix_reader=KmerFreqMatrixReader(),
        chain_type='IGH', strategy='NoKNeighbours'):
    assert chain_type in ['IGH', 'IGL', 'IGK']
    samples = [kmer_freq_matrix_reader.read(*x) for x in input_data]
    result = pd.concat((sample[strategy][chain_type] for sample in samples))
    return result


def concatenate_kmer_matrices_all_data(
        chain_type="IGH", strategy="NoKNeighbours",
        kmer_freq_matrix_reader=KmerFreqMatrixReader()):
    input_data = [("AbVitro/flu_time_course/FV/", ['25']),
                  ("AbVitro/flu_time_course/GMC/", ['8']),
                  ("AbVitro/flu_time_course/IDO/", []),
                  ("AbVitro/paired", [])]

    if chain_type == "IGH":
        input_data += [("age/", [])]

    result = concatenate_kmer_freq_matrices(
        input_data=input_data,
        kmer_freq_matrix_reader=kmer_freq_matrix_reader,
        chain_type=chain_type,
        strategy=strategy)
    return result
