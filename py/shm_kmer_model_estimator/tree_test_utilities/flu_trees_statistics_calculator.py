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
