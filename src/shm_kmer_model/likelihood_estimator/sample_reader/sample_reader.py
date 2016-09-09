import pandas as pd
import os


class SampleReader:
    """ Sample is a kmer frequency matrix.

    Reader tries the default folder structure:
    chain_type/indiv_number/strategy.

    Reader returns a dictionary of pandas.Panel.
    Items: numbers of datasets.
    Majo
    """
    def __init__(self):
        self.dir_kmer_statistics = 'kmer_statistics'
        self.filename_kmer_stats = 'kmer_statistics.csv'
        self.sep = ';'

    def read(self, root_dir, ignore_indiv_number=[],
             prefix_dir='/Users/andrewbzikadze/' +
                        'chihua/Sid/abzikadze/datasets/'):
        working_dir = os.path.join(prefix_dir, root_dir,
                                   self.dir_kmer_statistics)
        chain_types = os.listdir(working_dir)

        samples = {}
        for chain_type in chain_types:
            chain_type_dir = os.path.join(working_dir, chain_type)
            indiv_numbers = os.listdir(chain_type_dir)
            for indiv_number in indiv_numbers:
                if indiv_number in ignore_indiv_number:
                    continue
                indiv_number_dir = os.path.join(chain_type_dir, indiv_number)
                strategies = os.listdir(indiv_number_dir)
                for strategy in strategies:
                    strategy_path = os.path.join(indiv_number_dir, strategy)
                    if os.path.isdir(strategy_path):
                        sample_to_read = os.path.join(strategy_path,
                                                      self.filename_kmer_stats)
                        sample = pd.read_csv(sample_to_read,
                                             sep=self.sep,
                                             header=0,
                                             index_col=0)
                        try:
                            new_panel = \
                                pd.Panel({root_dir + indiv_number: sample})
                            samples[strategy][chain_type] = \
                                samples[strategy][chain_type].join(new_panel,
                                                                   'right')
                        except:
                            try:
                                samples[strategy][chain_type] = \
                                    pd.Panel({root_dir + indiv_number: sample})
                            except:
                                samples[strategy] = {chain_type: pd.Panel({
                                    root_dir + indiv_number: sample
                                })}
        return samples
