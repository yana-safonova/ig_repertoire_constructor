import os
import glob
from collections import defaultdict

from config.config import config
from chains.chains import Chains
from mutation_strategies.mutation_strategies import MutationStrategies


""" Sample is any data to be processed
that can be located by the default directory structure:
chain_type/indiv_number/strategy.
For example, a kmer frequency matrix.

Reader returns a dictionary of pandas.Panel.
Items: numbers of datasets.
"""
def read_samples(dir_data, filename_data, read_func,
                 prefix_dir, root_dir,
                 functionality,
                 ignore_indiv_number=[]):
    working_dir = os.path.join(prefix_dir, root_dir, dir_data)
    # default directory structure: chain_type/indiv_number/strategy
    pattern = os.path.join(working_dir, "*", "*", functionality, "*", filename_data)

    def rec_dd():
        return defaultdict(rec_dd)

    samples = rec_dd()
    for sample_filename in glob.iglob(pattern):
        sample_filename_spitted = sample_filename.split(os.sep)
        chain_type, indiv_number, _, strategy = sample_filename_spitted[-5:-1]
        chain_type, strategy = Chains[chain_type], MutationStrategies[strategy]
        if indiv_number in ignore_indiv_number:
            continue
        sample = read_func(sample_filename)
        samples[strategy][chain_type][indiv_number] = sample
    return samples
