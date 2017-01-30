import numpy as np

from kmer_utilities.kmer_utilities import central_nucl_indexes

class KmerMatrices(object):
    def __init__(self, fr_matrices, cdr_matrices):
        # self.fr_matrices = np.concatenate(fr_matrices.values(), 2)
        # self.cdr_matrices = np.concatenate(cdr_matrices.values(), 2)

        central_ind = central_nucl_indexes()
        central_ind = np.array(central_ind)
        matrix_shape = fr_matrices.itervalues().next().shape
        nonmut_matrix_shape = list(matrix_shape)
        nonmut_matrix_shape[1] -= 1
        central_ind_mask = np.zeros(matrix_shape, dtype=np.bool)

        central_ind_mask[xrange(len(central_ind_mask)),
                          central_ind] = True
        self.matrices = []
        for fr, cdr in zip(fr_matrices.values(), cdr_matrices.values()):
            fr_mut, cdr_mut = fr[central_ind_mask], cdr[central_ind_mask]
            fr_mut, cdr_mut = fr_mut[:, np.newaxis], cdr_mut[:, np.newaxis]
            fr_subst = fr[~central_ind_mask].reshape(nonmut_matrix_shape)
            cdr_subst = cdr[~central_ind_mask].reshape(nonmut_matrix_shape)
            subst = fr_subst + cdr_subst
            matrix = np.concatenate((fr_mut, cdr_mut, subst), 1)
            self.matrices.append(matrix[:, :, np.newaxis])
        self.matrices = np.concatenate(self.matrices, 2)


    def __str__(self):
        return self.matrices.__str__()

    def __repr__(self):
        return self.matrices.__repr__()
