import re
import numpy as np
import kmer_utilities.kmer_utilities as kmer_utilities


def hotspots():
    kmer_names = kmer_utilities.kmer_names()
    hotspot_pattern = r"^([AT][GA]C[ACT].)|(.[AGT]G[CT][AT])$"
    return filter(lambda kmer_name: re.search(hotspot_pattern, kmer_name),
                  kmer_names)


def coldspots():
    kmer_names = kmer_utilities.kmer_names()
    coldspot_pattern = r"^([CG][CT]C..)|(..G[GA][CG])$"
    return filter(lambda kmer_name: re.search(coldspot_pattern, kmer_name),
                  kmer_names)


def hotspots_indexes():
    kmer_names = kmer_utilities.kmer_names()
    return np.array([i for i in xrange(len(kmer_names))
                     if kmer_names[i] in hotspots()])


def coldspots_indexes():
    kmer_names = kmer_utilities.kmer_names()
    return np.array([i for i in xrange(len(kmer_names))
                     if kmer_names[i] in coldspots()])
