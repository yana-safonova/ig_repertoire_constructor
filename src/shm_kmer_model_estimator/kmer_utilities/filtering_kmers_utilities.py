#!/usr/bin/env python2

from __future__ import print_function

import kmer_utilities

import numpy as np
import pandas as pd


def filter_by_coverage(samples,
                       coverage_threshold=100,
                       mean_function=pd.DataFrame.median,
                       threshold_function=pd.Panel.max):
    """ Function filters a dataframe samples basing on the strategy
    samples -- normally a 3D array. 3rd dim -- for different samples.
    Reasonable choice for mean_function could be: np.median, np.mean
    ...                   threshold_function could be: np.min, np.max
    """
    coverage = mean_function(threshold_function(samples, axis=0), axis=1)
    ind_coverage = np.array(coverage > coverage_threshold)
    kmer_names = np.array(kmer_utilities.kmer_names())
    kmer_names = kmer_names[ind_coverage]
    return ind_coverage, samples.iloc[:, ind_coverage, :]
