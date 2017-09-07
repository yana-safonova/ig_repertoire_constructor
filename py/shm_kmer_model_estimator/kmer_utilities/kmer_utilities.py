import itertools
from memoize.memoized import memoized


def nucl_bases():
    return ['A', 'C', 'G', 'T']


@memoized
def kmer_names(kmer_len=5):
    """ Simple routine to get all kmer_names over DNA alphabet. """
    bases = nucl_bases()
    return [''.join(p) for p in itertools.product(bases, repeat=kmer_len)]


def central_nucl_indexes(kmer_len=5):
    """ Simple routine to get indexes in the alphabet
    of a central nucl in each kmer. """
    bases = nucl_bases()
    n_nucl = len(bases)
    half_kmer_len = kmer_len // 2
    indexes = [(x % n_nucl**(half_kmer_len + 1)) //
               n_nucl**half_kmer_len
               for x in xrange(len(nucl_bases())**kmer_len)]
    return indexes


@memoized
def nonmutated_index(kmer, kmer_len=5):
    """ Simple routine to get index in the alphabet
    of a central nucl in a kmer. """
    nonmutated_nucl = kmer[kmer_len // 2]
    return nucl_bases().index(nonmutated_nucl)


@memoized
def kmer_index(kmer, kmer_len=5):
    """ Simple routine to get index of a kmer in list of all kmers. """
    assert len(kmer) == kmer_len
    bases = nucl_bases()
    base = len(bases)
    p = 1
    index = 0
    for s in list(kmer)[::-1]:
        index += bases.index(s) * p
        p *= base
    return index


@memoized
def mutative_bases(kmer, kmer_len=5):
    """ Simple routine to get all bases to which central base can mutate """
    assert len(kmer) == kmer_len
    bases = nucl_bases()
    cent_n = kmer[kmer_len // 2]
    bases.remove(cent_n)
    return bases


@memoized
def reverse_mut_kmer(kmer, base, kmer_len=5):
    assert len(kmer) == kmer_len
    mut_kmer = list(kmer)
    mut_kmer[kmer_len // 2] = base
    return ''.join(mut_kmer)


def surround_kmers(i=None, kmer=None, bad_ind=[], kmer_len=5):
    if (i is None and kmer is None) or (i is not None and kmer is not None):
        raise ValueError("Supply either i or kmer and not both")
    if kmer is None:
        kmer = kmer_names()[i]
    surr_kmers, surr_kmer_inds = [], []
    bad_ind = set(bad_ind)
    for left in nucl_bases():
        for right in nucl_bases():
            surr_kmer = list(kmer)
            surr_kmer[0], surr_kmer[-1] = left, right
            surr_kmer = "".join(surr_kmer)
            surr_ind = kmer_index(surr_kmer)
            if surr_ind not in bad_ind:
                surr_kmer_inds.append(surr_ind)
                surr_kmers.append(surr_kmer)
    return surr_kmers, surr_kmer_inds
