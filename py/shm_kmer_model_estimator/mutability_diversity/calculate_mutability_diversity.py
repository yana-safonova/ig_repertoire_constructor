from __future__ import division
import itertools

import numpy as np
import pandas as pd


def calculate_mutability_full(kmer_matrices):
    matrices = kmer_matrices.matrices
    full = np.sum(matrices[:, :4, :], axis = 1)
    mut = matrices[:, 0, :] + matrices[:, 2, :]
    return mut / full


def calculate_mutability_fr(kmer_matrices):
    matrices = kmer_matrices.matrices
    return matrices[:, 0, :] / (matrices[:, 0, :] + matrices[:, 1, :])


def calculate_mutability_cdr(kmer_matrices):
    matrices = kmer_matrices.matrices
    return matrices[:, 2, :] / (matrices[:, 2, :] + matrices[:, 3, :])


def calculate_substitution(kmer_matrices):
    matrices = kmer_matrices.matrices
    norm = np.sum(matrices[:, 4:, :], axis = 1)
    norm = norm[:, np.newaxis, :]
    return matrices[:, 4:, :] / norm
