import numpy as np

from kmer_utilities.kmer_utilities import central_nucl_indexes, \
                                          kmer_names, \
                                          kmer_index

class KmerMatrices(object):
    def __init__(self, fr_matrices=None, cdr_matrices=None):
        if fr_matrices is None or cdr_matrices is None:
            return

        matrix_shape = fr_matrices.itervalues().next().shape
        nonmut_matrix_shape = list(matrix_shape)
        nonmut_matrix_shape[1] -= 1

        central_ind = central_nucl_indexes()
        central_ind = np.array(central_ind)
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

        self.fr_matrices = \
            [fr[:, :, np.newaxis] for fr in fr_matrices.itervalues()]
        self.cdr_matrices = \
            [cdr[:, :, np.newaxis] for cdr in cdr_matrices.itervalues()]
        self.fr_matrices = np.concatenate(self.fr_matrices, 2)
        self.cdr_matrices = np.concatenate(self.cdr_matrices, 2)

    @classmethod
    def FromKmerMatricesList(cls, kmer_matrices):
        obj = cls()
        obj.matrices = \
            np.concatenate([x.matrices for x in kmer_matrices], 2)
        obj.fr_matrices = \
            np.concatenate([x.fr_matrices for x in kmer_matrices], 2)
        obj.cdr_matrices = \
            np.concatenate([x.cdr_matrices for x in kmer_matrices], 2)
        return obj

    def __str__(self):
        return self.matrices.__str__()

    def __repr__(self):
        return self.matrices.__repr__()

    def __getitem__(self, key):
        if isinstance(key, str):
            if key not in kmer_names():
                raise TypeError, "Not a kmer string argument: %s" % key
            index = kmer_index(key)
            return self.matrices[index, :, :].T
        if isinstance(key, int):
            return self.matrices[:, :, key]
        if isinstance(key, slice):
            return self.matrices[:, :, key.start:key.stop:key.step]
        raise TypeError, "Unknown key type: %s" % key

