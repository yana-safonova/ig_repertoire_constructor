#!/usr/bin/env python2

from __future__ import print_function

import numpy as np
import pandas as pd

import mismatch_finder.mismatch_finder as mismatch_finder

class LikehoodCalculator(object):
    def __init__(self, model, kmer_len=5,
                 mismatch_finder=mismatch_finder.NoKNeighboursMismatchFinder(5)):
        self.model = model
        self.kmer_len = kmer_len
        self.half_kmer_len = self.kmer_len // 2
        self.mismatch_finder = mismatch_finder

    def calculate_likelihood(self, source, destination):
        assert len(source) == len(destination)
        mismatch_positions = self.mismatch_finder.find_mismatch_positions(source, destination)
        likelihood = 0
        for mismatch_position in mismatch_positions:
            source_fivemer = source[(mismatch_position - self.half_kmer_len) : \
                                    (mismatch_position + self.half_kmer_len + 1)]
            likelihood += np.log(self.model.expect_mut_prob(source_fivemer)) + \
                          np.log(self.model.expect_subst_prob(source_fivemer,
                                                              destination[mismatch_position]))
        return likelihood
