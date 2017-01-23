#!/usr/bin/env python2

from __future__ import print_function

import numpy as np
from abc import ABCMeta, abstractmethod


class AbstractMismatchFinder:
    __metaclass__ = ABCMeta

    @abstractmethod
    def find_mismatch_positions(self, seq1, seq2):
        pass


class TrivialMismatchFinder(AbstractMismatchFinder):
    def __init__(self, config):
        pass

    def find_mismatch_positions(self, seq1, seq2):
        if len(seq1) != len(seq2):
            raise ValueError('Lengths of seq1 and seq2')
        return [i for i in xrange(len(seq1)) if seq1[i] != seq2[i]]


class NoKNeighboursMismatchFinder(AbstractMismatchFinder):
    def __init__(self, config):
        self.k_mer_len = config.k_mer_len

    def find_mismatch_positions(self, seq1, seq2):
        assert len(seq1) == len(seq2)
        semi_len = self.k_mer_len // 2

        mask = np.array([a != b for a, b in zip(seq1, seq2)])
        conv = np.convolve(mask, np.ones(self.k_mer_len))
        conv = conv[semi_len:-semi_len]

        mask[:semi_len] = 0
        mask[-semi_len:] = 0
        return np.where((conv == 1) & (mask))[0].tolist()
