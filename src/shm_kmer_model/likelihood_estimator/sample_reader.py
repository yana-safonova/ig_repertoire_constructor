from pandas import read_csv
import numpy as np
import os

class SampleReader:
    def __init__(self):
        self.dir_kmer_statistics = 'kmer_statistics'
        self.filename_kmer_statistics = 'kmer_statistics.csv'
        self.sep = ';'
    
    def read(self, root_dir, ignore_indiv_number=[]):
        working_dir = os.path.join(root_dir, self.dir_kmer_statistics)
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
                                                      self.filename_kmer_statistics)
                        sample = read_csv(sample_to_read,
                                          sep=self.sep,
                                          header=0,
                                          index_col=0)
                        try:
                            samples[strategy][chain_type] = \
                                np.concatenate((samples[strategy][chain_type], np.array(sample, ndmin=3)))
                        except:
                            try:
                                samples[strategy][chain_type] = np.array(sample, ndmin=3)
                            except:
                                samples[strategy] = {chain_type: np.array(sample, ndmin=3)}
        return samples
