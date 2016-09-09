import itertools

import numpy as np
import pandas as pd


def calculate_mutability_diversity(samples,
                                   data_to_plot,
                                   good_coverage_indices=np.arange(4**5),
                                   kmer_names=np.array([''.join(p) for p in
                                                        itertools.product(
                                                        ['A', 'C', 'G', 'T'],
                                                        repeat=5)]),
                                   nonmutated_ind=(np.arange(4**5) / 4**2)
                                   % 4):
    samples = np.copy(samples)
    good_coverage_indices = np.copy(good_coverage_indices)

    kmer_names = kmer_names[good_coverage_indices]
    nonmutated_ind = nonmutated_ind[good_coverage_indices]

    order_array = np.array([np.median(data_to_plot(samples[:, i, :], i,
                                                   nonmutated_ind))
                            for i in xrange(samples.shape[1])])
    ind = np.argsort(-order_array)

    nonmutated_ind = nonmutated_ind[ind]
    samples = samples[:, ind, :]

    if data_to_plot == calculate_full_substitution:
        mutability_dataframe = pd.Panel(
            np.array([data_to_plot(samples[:, i, :], i, nonmutated_ind)
                      for i in xrange(samples.shape[1])]),
            items=kmer_names[ind])
    else:
        mutability_dataframe = pd.DataFrame(
            np.array([data_to_plot(samples[:, i, :], i, nonmutated_ind)
                      for i in xrange(samples.shape[1])]),
            index=kmer_names[ind])
        mutability_dataframe = mutability_dataframe.T

    return ind, mutability_dataframe


def calculate_mutability(row, row_index, nonmutated_ind):
    return 1 - row[:, nonmutated_ind[row_index]] / np.sum(row, 1, dtype=float)


def calculate_substitution(row, row_index, nonmutated_ind, i):
    return row[:, (nonmutated_ind[row_index] + i) % 4] / \
           (np.sum(row, 1, dtype=float) - row[:, nonmutated_ind[row_index]])


def calculate_substitution_1(row, row_index, nonmutated_ind):
    return calculate_substitution(row, row_index, nonmutated_ind, 1)


def calculate_substitution_2(row, row_index, nonmutated_ind):
    return calculate_substitution(row, row_index, nonmutated_ind, 2)


def calculate_substitution_3(row, row_index, nonmutated_ind):
    return calculate_substitution(row, row_index, nonmutated_ind, 3)


def calculate_full_substitution(row, row_index, nonmutated_ind):
    return [calculate_substitution(row, row_index, nonmutated_ind, i)
            for i in xrange(3)]
