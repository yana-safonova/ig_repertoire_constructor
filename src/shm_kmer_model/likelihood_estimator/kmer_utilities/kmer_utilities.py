import itertools

def nucl_bases():
    return ['A', 'C', 'G', 'T']

def kmer_names(kmer_len=5):
    bases = nucl_bases()
    return [''.join(p) for p in itertools.product(bases, repeat=kmer_len)]

def central_nucl_indexes(kmer_len=5, debug=False):
    bases = nucl_bases()
    n_nucl = len(bases)
    half_kmer_len = kmer_len // 2
    indexes = [(x % n_nucl**(half_kmer_len + 1)) //\
               n_nucl**half_kmer_len \
               for x in xrange(len(nucl_bases())**kmer_len)]
    if debug:
        names = kmer_names()
        for i in xrange(len(bases)**kmer_len):
            assert names[i][kmer_len // 2] == bases[indexes[i]]
    return indexes
