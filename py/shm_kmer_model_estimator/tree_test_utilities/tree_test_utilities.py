#!/usr/bin/env python2

from __future__ import print_function

import numpy as np
import pandas as pd

import likelihood_calculator.likelihood_calculator as likelihood_calculator
import mismatch_finder.mismatch_finder as mismatch_finder

from disk_memoize.disk_memoize import memoize_to_disk


class TreeTestResults(object):
    def get_accuracy__(self, res):
        return np.nanmean([x[0] < x[1] for x in res])

    def __init__(self, lklhs):
        self.lklhs = lklhs
        self.accuracies = np.array([self.get_accuracy__(lklh) for lklh in lklhs])
        self.full_accuracy = self.get_accuracy__(np.concatenate(lklhs)) if len(lklhs) else np.nan


class TreeTester(object):
    def __init__(self, minimal_size_filtered_tree=20):
        self.minimal_size_filtered_tree = minimal_size_filtered_tree

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
        for source, destination in \
                tree[['Src_CDR3', 'Dst_CDR3']].itertuples(False):
            numb_mis_pos = len(mf.find_mismatch_positions(source, destination))
            mutated_indexes.append(numb_mis_pos != 0)
        tree = tree.loc[mutated_indexes]
        return tree

    def get_trees(self, tree_paths, mismatch_strategy='Trivial'):
        trees = [self.__read_tree(tree_path) for tree_path in tree_paths]
        trees = [self.__filter_tree(tree, mismatch_strategy) for tree in trees]
        trees = [tree for tree in trees if len(tree)]
        return trees

    def get_likelihood_statistics_tree(self, model, tree_path=None,
                                       mismatch_strategy='Trivial'):
        tree = self.__read_tree(tree_path)

        tree = self.__filter_tree(tree, mismatch_strategy)
        if tree.shape[0] < self.minimal_size_filtered_tree:
            return np.array([], dtype=np.dtype('float, float'))

        lkhd_calc = likelihood_calculator.LikelihoodCalculator(model)

        results = []
        for source, dest in \
                tree[['Src_CDR3', 'Dst_CDR3']].itertuples(False):
            results.append((lkhd_calc.calculate_loglikelihood(source, dest),
                            lkhd_calc.calculate_loglikelihood(dest, source)))
        return np.array(results, dtype=np.dtype('float, float'))

    def get_likelihood_statistics_trees(self, model, tree_paths,
                                        mismatch_strategy='Trivial'):
        results = [self.get_likelihood_statistics_tree(model, tree_path, mismatch_strategy)
                   for tree_path in tree_paths]
        results = [x for x in results if x.size]
        return TreeTestResults(results)
