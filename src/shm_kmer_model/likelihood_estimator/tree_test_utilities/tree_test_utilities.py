#!/usr/bin/env python2

from __future__ import print_function

import numpy as np
import pandas as pd

import likelihood_calculator.likelihood_calculator as likelihood_calculator
import mismatch_finder.mismatch_finder as mismatch_finder

class TreeTester(object):
    def __init__(self, minimal_size_filtered_tree=20):
        self.minimal_size_filtered_tree=minimal_size_filtered_tree

    def __read_tree(self, tree_path):
        tree = pd.read_csv(tree_path, sep='\t')
        return tree

    def __filter_tree(self, tree, mismatch_strategy):
        tree = tree.copy()
        tree = tree.loc[(tree.Edge_type == 'directed')]
        tree.set_index(np.arange(tree.shape[0]), inplace=True)

        if mismatch_strategy == 'Trivial':
            mf = mismatch_finder.TrivialMismatchFinder()
        elif mismatch_strategy == 'NoKNeighbours':
            mf = mismatch_finder.NoKNeighboursMismatchFinder(5)
        else:
            raise ValueError('Not supported mismatch_strategy')

        mutated_indexes = []
        for source, destination in tree[['Src_CDR3', 'Dst_CDR3']].itertuples(False):
            mutated_indexes.append(len(mf.find_mismatch_positions(source, destination)) != 0)
        tree = tree.loc[mutated_indexes]
        return tree

    def get_consistency_statistics(self, model, tree_path=None, tree=None,
                                   mismatch_strategy='NoKNeighbours'):
        if tree_path is None and tree is None:
            raise ValueError('tree itself or tree_path must be supplied')
        if tree is None:
            tree = self.__read_tree(tree_path)

        tree = self.__filter_tree(tree, mismatch_strategy)
        if tree.shape[0] < self.minimal_size_filtered_tree:
            return np.array([])

        lkhd_calc = likelihood_calculator.LikehoodCalculator(model)

        results = []
        for source, destination in tree[['Src_CDR3', 'Dst_CDR3']].itertuples(False):
            results += [lkhd_calc.calculate_likelihood(source, destination) < \
                        lkhd_calc.calculate_likelihood(destination, source)]

        return np.array(results)
