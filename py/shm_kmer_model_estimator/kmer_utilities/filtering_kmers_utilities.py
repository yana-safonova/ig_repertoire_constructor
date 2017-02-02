#!/usr/bin/env python2

from __future__ import print_function
import numpy as np

import kmer_utilities
from kmer_matrices.kmer_matrices import KmerMatrices


def filter_by_coverage(kmer_matrices,
                       coverage_threshold=100,
                       mean_function=np.median,
                       threshold_function=np.max):
    """ Function filters kmer_matrices basing on the strategy
    kmer_matrices -- a KmerMatrices class object

    mean_function -- across samples.
    threshold_function -- across nucleotides.

    Reasonable choice for mean_function could be: np.median, np.mean
    ...                   threshold_function could be: np.min, np.max
    """
    pure_matrices = kmer_matrices.matrices
    coverage = mean_function(threshold_function(pure_matrices, axis=1), axis=1)
    covered_ind  = np.array(coverage > coverage_threshold)
    return KmerMatrices.FromKmerMatricesByFiltering(kmer_matrices, covered_ind)
