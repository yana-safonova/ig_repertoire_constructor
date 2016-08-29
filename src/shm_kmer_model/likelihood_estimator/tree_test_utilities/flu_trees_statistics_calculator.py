#!/usr/bin/env python2

from __future__ import print_function

import numpy as np
import pandas as pd

from special_utils.cd_manager import cd
import shm_kmer_model.shm_kmer_model as shm_kmer_model
import shm_kmer_model_estimator.standard_model_estimations as standard_model_estimations


class FluTreesStatisticsCalculator(object):
    def __init__(self, n_largest_trees=5,
                 data_path_prefix=\
                 '/Users/andrewbzikadze/chihua/home/aslabodkin/trees/trees_res',
                 yale_model=shm_kmer_model.YaleSHM_Model(),
                 cab_model=standard_model_estimations.full_igh_model()):
        self.n_largest_trees = n_largest_trees
        self.data_path_prefix = data_path_prefix
        self.yale_model = yale_model
        self.cab_model = cab_model
        self.flu_ind_names = ['IDO', 'FV', 'GMC']
        self.strategies = ['Trivial', 'NoKNeighbours']
        self.model_names = ['Yale', 'CAB_NoKNeighbours', 'CAB_Trivial']

    def get_flu_trees_paths(self):
        from special_utils.largest_files_in_dir import n_largest_files
        import os, re
        paths = []
        for flu_ind in self.flu_ind_names:
            ind_prefix = os.path.join(self.data_path_prefix, flu_ind)
            with cd(ind_prefix):
                dir_list = filter(os.path.isdir, os.listdir(ind_prefix))
                dir_list = filter(lambda x: re.match(r'.*_heavy', x), dir_list)
                for dir_id in dir_list:
                    dataset_prefix = os.path.join(os.getcwd(), dir_id, 'clonal_trees')
                    paths += n_largest_files(dataset_prefix, self.n_largest_trees)
        return paths

    def get_flu_likelihood_statistics(self, tester, model_mode):
        results = dict.fromkeys(self.model_names)
        for key in results:
            results[key] = {self.strategies[0]: [], self.strategies[1]: []}

        def append_lkhd_stat(dataset, model):
            dataset.append(tester.get_likelihood_statistics(model,
                mismatch_strategy=strategy,
                tree_path=flu_tree_path,
                model_mode=model_mode))
            return dataset

        flu_trees_paths = self.get_flu_trees_paths()
        for flu_tree_path in flu_trees_paths:
            for strategy in self.strategies:
                results['Yale'][strategy] = \
                    append_lkhd_stat(results['Yale'][strategy],
                        self.yale_model)
                results['CAB_NoKNeighbours'][strategy] = \
                    append_lkhd_stat(results['CAB_NoKNeighbours'][strategy],
                        self.cab_model['NoKNeighbours']['IGH'])
                results['CAB_Trivial'][strategy] = \
                    append_lkhd_stat(results['CAB_Trivial'][strategy],
                        self.cab_model['Trivial']['IGH'])

        for k1 in results:
            for k2 in results[k1]:
                results[k1][k2] = np.array(results[k1][k2])

        return results
