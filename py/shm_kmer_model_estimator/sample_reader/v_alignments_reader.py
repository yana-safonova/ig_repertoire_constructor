import os

import special_utils.os_utils as os_utils
import v_alignments_reader_utilities.alignment_reader as alignment_reader


class VAlignmentsReader(object):
    def __init__(self,
                 dir_data='vjf_output_on_source',
                 filename_data='valignments.fa'):
        self.dir_data = dir_data
        self.filename_data = filename_data
        self.reader = alignment_reader.GermlineAlignmentReader()
        self.read_func = (lambda x: self.reader.read(x))

    def read(self, root_dir, ignore_indiv_number=[],
             prefix_dir='/Users/andrewbzikadze/' +
                        'chihua/Sid/abzikadze/datasets/'):
        working_dir = os.path.join(prefix_dir, root_dir,
                                   self.dir_data)
        chain_types = os_utils.list_only_dirs(working_dir)

        read_matrices = {}
        for chain_type in chain_types:
            chain_type_dir = os.path.join(working_dir, chain_type)
            indiv_numbers = os_utils.list_only_dirs(chain_type_dir)
            for indiv_number in indiv_numbers:
                if indiv_number in ignore_indiv_number:
                    continue
                indiv_number_dir = os.path.join(chain_type_dir, indiv_number)
                print(indiv_number_dir)
                if os.path.isdir(indiv_number_dir):
                    matrix_to_read = os.path.join(indiv_number_dir,
                                                  self.filename_data)
                    matrix = self.read_func(matrix_to_read)
                    new_key = root_dir + indiv_number_dir
                    try:
                        read_matrices[chain_type][new_key] = matrix
                    except:
                        read_matrices[chain_type] = {new_key: matrix}
        return read_matrices
