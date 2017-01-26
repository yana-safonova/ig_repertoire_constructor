""" This module contains standard samples from AbVitro and age datasets.
Func concatenate_samples creates a Panel that is as full as possible. """

from kmer_freq_matrix_reader import KmerFreqMatrixReader

import pandas as pd

kmer_freq_matrix_reader = KmerFreqMatrixReader()

kmer_freq_matrices_fv     = kmer_freq_matrix_reader.read('AbVitro/flu_time_course/FV/', ['25'])
kmer_freq_matrices_gmc    = kmer_freq_matrix_reader.read('AbVitro/flu_time_course/GMC/', ['8'])
kmer_freq_matrices_ido    = kmer_freq_matrix_reader.read('AbVitro/flu_time_course/IDO/')
kmer_freq_matrices_age    = kmer_freq_matrix_reader.read('age/')
kmer_freq_matrices_paired = kmer_freq_matrix_reader.read('AbVitro/paired/')


def concatenate_kmer_freq_matrices(chain_type='IGH', strategy='NoKNeighbours'):
    assert chain_type in ['IGH', 'IGL', 'IGK']
    result = pd.concat((kmer_freq_matrices_fv[strategy][chain_type],
                        kmer_freq_matrices_ido[strategy][chain_type],
                        kmer_freq_matrices_gmc[strategy][chain_type],
                        kmer_freq_matrices_paired[strategy][chain_type]))
    if chain_type != 'IGH':
        return result
    else:
        return pd.concat((result, kmer_freq_matrices_age[strategy][chain_type]))
