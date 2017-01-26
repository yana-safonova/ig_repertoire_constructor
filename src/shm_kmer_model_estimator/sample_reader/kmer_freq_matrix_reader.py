import pandas as pd
from sample_reader import SampleReader
from disk_memoize.disk_memoize import memoize_to_disk


# @memoize_to_disk
def read_matrix(x):
    return pd.read_csv(x, sep=';', header=0, index_col=0)


class KmerFreqMatrixReader(SampleReader):
    def __init__(self, dir_data='kmer_statistics',
                 filename_data='kmer_statistics.csv',
                 read_func=read_matrix):
        super(self.__class__, self).__init__(dir_data,
                                             filename_data,
                                             read_func)

    def __post_process(self, read_matrix):
        for strategy in read_matrix:
            for chain_type in read_matrix[strategy]:
                read_matrix[strategy][chain_type] = \
                    pd.Panel.from_dict(read_matrix[strategy][chain_type],
                                       intersect=True)
        return read_matrix

    def read(self, root_dir, ignore_indiv_number=[],
             prefix_dir='/Users/andrewbzikadze/' +
                        'chihua/Sid/abzikadze/datasets/'):
        read_matrices = super(self.__class__, self).read(root_dir,
                                                         ignore_indiv_number,
                                                         prefix_dir)
        read_matrices = self.__post_process(read_matrices)
        return read_matrices
