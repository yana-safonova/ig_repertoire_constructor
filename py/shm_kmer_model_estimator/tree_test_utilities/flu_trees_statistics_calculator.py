#!/usr/bin/env python2

from __future__ import print_function

import os
import re

import numpy as np

from special_utils.cd_manager import cd
from special_utils.largest_files_in_dir import n_largest_files
from shm_kmer_model import cab_shm_model, yale_model
import shm_kmer_model_estimator.standard_model_estimations as \
    standard_model_estimations
from chains import chains

from disk_memoize.disk_memoize import memoize_to_disk 


def get_flu_trees_paths(chain_type=chains.Chains.IGH, n_largest_trees=5,
                         data_path_prefix='/home/aslabodkin/antevolo_project/trees_res'):
    chain_path = {chains.Chains.IGH: 'heavy',
                  chains.Chains.IGL: 'lambda',
                  chains.Chains.IGK: 'kappa'}[chain_type]

    return get_flu_trees_paths__(chain_path, n_largest_trees, data_path_prefix)


@memoize_to_disk
def get_flu_trees_paths__(chain_path, n_largest_trees, data_path_prefix):
    paths = []
    flu_ind_names = ['IDO', 'FV', 'GMC']
    for flu_ind in flu_ind_names:
        ind_prefix = os.path.join(data_path_prefix, flu_ind)
        with cd(ind_prefix):
            dir_list = filter(os.path.isdir, os.listdir(ind_prefix))
            dir_list = filter(lambda x: re.match(r'.*_' + chain_path, x),
                              dir_list)
            for dir_id in dir_list:
                dataset_prefix = os.path.join(os.getcwd(), dir_id,
                                              'clonal_trees')
                paths += n_largest_files(dataset_prefix,
                                         n_largest_trees)
    return paths


class FluTreesStatisticsCalculator(object):
    def __init__(self, strategy, chain, n_largest_trees=5,
                 data_path_prefix='/Users/andrewbzikadze/chihua/home/' +
                 'aslabodkin/antevolo_project/trees_res'):
                 #full_igh_model_extended_base()):
        self.n_largest_trees = n_largest_trees
        self.data_path_prefix = data_path_prefix
        self.yale_model = yale_model.YaleSHM_Model
        self.cab_model = cab_shm_model.CAB_SHM_Model(strategy, chain)
        self.flu_ind_names = ['IDO', 'FV', 'GMC']
        self.model_names = ['Yale', 'CAB_NoKNeighbours', 'CAB_Trivial']

    def get_flu_likelihood_statistics(self, tester, model_mode,
                                      chain_type='IGH'):
        results = dict.fromkeys(self.model_names)
        for key in results:
            results[key] = {self.strategies[0]: [], self.strategies[1]: []}

        def append_lkhd_stat(dataset, model):
            dataset.append(tester.get_likelihood_statistics(model,
                           mismatch_strategy=strategy,
                           tree_path=flu_tree_path,
                           model_mode=model_mode))
            return dataset

        flu_trees_paths = self.get_flu_trees_paths(chain_type)
        for flu_tree_path in flu_trees_paths:
            for strategy in self.strategies:
                results['Yale'][strategy] = \
                    append_lkhd_stat(results['Yale'][strategy],
                                     self.yale_model)
                results['CAB_NoKNeighbours'][strategy] = \
                    append_lkhd_stat(results['CAB_NoKNeighbours'][strategy],
                                     self.cab_model['NoKNeighbours'][chain_type])
                results['CAB_Trivial'][strategy] = \
                    append_lkhd_stat(results['CAB_Trivial'][strategy],
                                     self.cab_model['Trivial'][chain_type])
                results['CAB_NoKNeighbours_extended_base'][strategy] = \
                    append_lkhd_stat(results['CAB_NoKNeighbours_extended_base'][strategy],
                                     self.cab_model_extended_base['NoKNeighbours'][chain_type])
                results['CAB_Trivial_extended_base'][strategy] = \
                    append_lkhd_stat(results['CAB_Trivial_extended_base'][strategy],
                                     self.cab_model_extended_base['Trivial'][chain_type])

        for k1 in results:
            for k2 in results[k1]:
                results[k1][k2] = np.array(results[k1][k2])

        return results
