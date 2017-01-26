import os

import special_utils.os_utils as os_utils


class SampleReader(object):
    """ Sample is any data to be processed
    that can be located by the default folder structure:
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
        chain_types = os_utils.list_only_dirs(working_dir)

        samples = {}
        for chain_type in chain_types:
            chain_type_dir = os.path.join(working_dir, chain_type)
            indiv_numbers = os_utils.list_only_dirs(chain_type_dir)
            for indiv_number in indiv_numbers:
                if indiv_number in ignore_indiv_number:
                    continue
                indiv_number_dir = os.path.join(chain_type_dir, indiv_number)
                strategies = os_utils.list_only_dirs(indiv_number_dir)
                for strategy in strategies:
                    strategy_path = os.path.join(indiv_number_dir, strategy)
                    if os.path.isdir(strategy_path):
                        sample_to_read = os.path.join(strategy_path,
                                                      self.filename_data)
                        sample = self.read_func(sample_to_read)
                        new_key = root_dir + indiv_number_dir
                        try:
                            samples[strategy][chain_type][new_key] = \
                                    sample
                        except:
                            try:
                                samples[strategy][chain_type] = \
                                    {new_key: sample}
                            except:
                                samples[strategy] = {chain_type: {new_key:
                                                                  sample}}
        return samples
