#!/usr/bin/env python2

from __future__ import print_function

import numpy as np
from abc import ABCMeta, abstractmethod


class AbstractMismatchFinder:
    """ An abstract class for mismatch finding strategy.
    Two subclasses are currently available, basing on stragegies:
    trivial and nokneighbors. """
    __metaclass__ = ABCMeta

    def __init__(self, kmer_len=5):
        self.kmer_len = kmer_len
        self.half_kmer_len = self.kmer_len // 2

    @abstractmethod
    def find_mismatch_positions(self, seq1, seq2):
        pass

    def filter_N_in_mismatch_positions(self, seq1, seq2,
            mismatch_positions):
        extract_kmer = (lambda seq, pos: seq[(pos - self.half_kmer_len):\
                                             (pos + self.half_kmer_len + 1)])
        filter_condition = (lambda pos: 'N' != seq2[pos] and \
                                        'N' not in extract_kmer(seq1, pos))
        return filter(filter_condition, mismatch_positions)


class TrivialMismatchFinder(AbstractMismatchFinder):
    """ Trivial strategy considers __all__ mismatches. """
    def __init__(self, kmer_len=5):
        super(TrivialMismatchFinder, self).__init__(kmer_len)

    def find_mismatch_positions(self, seq1, seq2):
        if len(seq1) != len(seq2):
            raise ValueError('Lengths of seq1 and seq2')
        mismatches_range = xrange(self.half_kmer_len,
                                  len(seq1) - self.half_kmer_len)
        mismatch_positions = [i for i in mismatches_range if seq1[i] != seq2[i]]
        mismatch_positions = \
                super(TrivialMismatchFinder, self).\
                filter_N_in_mismatch_positions(seq1, seq2, mismatch_positions)
        return mismatch_positions


class NoKNeighboursMismatchFinder(AbstractMismatchFinder):
    """ NoKNeighbours strategy considers only mutations that are
    far enough from each other. """
    def __init__(self, kmer_len=5):
        super(NoKNeighboursMismatchFinder, self).__init__(kmer_len)

    def find_mismatch_positions(self, seq1, seq2):
        if len(seq1) != len(seq2):
            raise ValueError('Lengths of seq1 and seq2')
        semi_len = self.kmer_len // 2

        mask = np.array([a != b for a, b in zip(seq1, seq2)])
        conv = np.convolve(mask, np.ones(self.kmer_len))
        conv = conv[semi_len:-semi_len]

        mask[:semi_len] = 0
        mask[-semi_len:] = 0
        mismatch_positions = np.where((conv == 1) & (mask))[0].tolist()
        mismatch_positions = \
                super(NoKNeighboursMismatchFinder, self).\
                filter_N_in_mismatch_positions(seq1, seq2, mismatch_positions)
        return mismatch_positions
