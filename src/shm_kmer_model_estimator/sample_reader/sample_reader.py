import os
import glob
from collections import defaultdict


class SampleReader(object):
    """ Sample is any data to be processed
    that can be located by the default directory structure:
    chain_type/indiv_number/strategy.
    For example, a kmer frequency matrix.

    Reader returns a dictionary of pandas.Panel.
    Items: numbers of datasets.
    """
    def __init__(self, dir_data, filename_data, read_func):
        self.dir_data = dir_data
        self.filename_data = filename_data
        self.read_func = read_func

    def read(self, root_dir, ignore_indiv_number=[],
             prefix_dir='/Users/andrewbzikadze/' +
                        'chihua/Sid/abzikadze/datasets/'):
        working_dir = os.path.join(prefix_dir, root_dir,
                                   self.dir_data)
        # default directory structure: chain_type/indiv_number/strategy
        pattern = os.path.join(working_dir, "*", "*", "*", self.filename_data)

        def rec_dd():
            return defaultdict(rec_dd)

        samples = rec_dd()
        for sample_filename in glob.iglob(pattern):
            sample_filename_spitted = sample_filename.split(os.sep)
            chain_type, indiv_number, strategy = sample_filename_spitted[-4:-1]
            if indiv_number in ignore_indiv_number:
                continue
            sample = self.read_func(sample_filename)
            samples[strategy][chain_type][indiv_number] = sample
        return samples
