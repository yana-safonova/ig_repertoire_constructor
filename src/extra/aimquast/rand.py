


def comb2(n):
    return n * (n - 1) // 2


def randstats(X, Y):
    from collections import defaultdict

    m = defaultdict(int)
    mX = defaultdict(int)
    mY = defaultdict(int)

    for x, y in zip(X, Y):
        m[(x, y)] += 1
        mX[x] += 1
        mY[y] += 1

    S00 = S01 = S10 = S11 = 0

    assert len(Y) == len(X)
    N = len(X)

    for x, y in zip(X, Y):
        _m = m[(x, y)]
        _mX = mX[x]
        _mY = mY[y]
        S00 += _m - 1
        S01 += _mX - _m
        S10 += _mY - _m
        S11 += N - _mX - _mY + _m


    S00 //= 2
    S01 //= 2
    S10 //= 2
    S11 //= 2

    assert S00 + S11 + S01 + S10 == N * (N - 1) // 2

    return S00, S11, S01, S10


def rand_adj(X, Y):
    from collections import defaultdict

    m = defaultdict(int)
    mX = defaultdict(int)
    mY = defaultdict(int)

    for x, y in zip(X, Y):
        m[(x, y)] += 1
        mX[x] += 1
        mY[y] += 1

    S00 = S01 = S10 = S11 = 0

    assert len(Y) == len(X)
    N = len(X)

    for x, y in zip(X, Y):
        _m = m[(x, y)]
        _mX = mX[x]
        _mY = mY[y]
        S00 += _m - 1
        S01 += _mX - _m
        S10 += _mY - _m
        S11 += N - _mX - _mY + _m

    S00 //= 2
    S01 //= 2
    S10 //= 2
    S11 //= 2


    # Compute the ARI using the contingency data
    sum_comb_c = sum(comb2(_) for _ in mY.itervalues())
    sum_comb_k = sum(comb2(_) for _ in mX.itervalues())
    sum_comb = sum(comb2(_) for _ in m.itervalues())

    assert sum_comb == S00
    assert S01 == sum_comb_k - S00
    assert S10 == sum_comb_c - S00
    assert S11 == comb2(N) - S01 - S10 - S00


    prod_comb = (sum_comb_c * sum_comb_k) / float(comb2(N))
    mean_comb = (sum_comb_k + sum_comb_c) / 2.

    return (sum_comb - prod_comb) / (mean_comb - prod_comb)


def purity(X, Y):
    # TODO
    pass


def normalized_mutual_information(X, Y): #  NMI
    pass


def rand(X, Y):
    S00, S11, S01, S10 = randstats(X, Y)

    return float(S00 + S11) / float(S00 + S11 + S10 + S01)

def FM_index(X, Y):
    import math
    S00, S11, S01, S10 = randstats(X, Y)

    return float(S00) / math.sqrt((S00 + S10) * (S00 + S01))

def jaccard(X, Y):
    S00, S11, S01, S10 = randstats(X, Y)

    return float(S00) / (S00 + S10 + S01)

assert rand([0, 1], [0, 1]) == 1

import sklearn.metrics



assert sklearn.metrics.adjusted_rand_score([0, 1, 1], [0, 1, 2]) ==  rand_adj([0, 1, 1], [0, 1, 2])
assert sklearn.metrics.adjusted_rand_score([0, 1, 1], [0, 1, 1]) == rand_adj([0, 1, 1], [0, 1, 1])
assert sklearn.metrics.adjusted_rand_score([0, 0, 1, 2], [0, 0, 1, 1]) == rand_adj([0, 0, 1, 2], [0, 0, 1, 1])




def parse_rcm(filename):
    rcm = {}
    with open(filename) as f:
        for line in f:
            id, cluster = line.strip().split("\t")
            id = id.strip()
            cluster = cluster.strip()
            rcm[id] = cluster

    return rcm

def rcm_vs_rcm(rcm1, rcm2, score=FM_index):
    rcm1, rcm2 = parse_rcm(rcm1), parse_rcm(rcm2)

    clustering1 = []
    clustering2 = []
    for id, cluster1 in rcm1.iteritems():
        if id in rcm2:
            cluster2 = rcm2[id]
            clustering1.append(cluster1)
            clustering2.append(cluster2)

    return score(clustering1, clustering2)

def make_ideal_rcm(filename, output):
    import re

    with open(filename) as fh, open(output, "w") as out:
        for line in fh:
            id, cluster = line.strip().split("\t")
            cluster = int(cluster)
            m = re.match("^antibody_(\\d+)_", id)
            ant = m.groups()[0]

            out.write("%s\t%s\n" % (id, ant))


def parse_cluster_mult(id):
    import re
    id = str(id)
    m = re.match(r"^cluster___([\dA-Za-z]+)___size___(\d+)$", id)
    if m:
        g = m.groups()
        cluster = g[0]
        mult = int(g[1])
        return cluster, mult
    else:
        return None

def idFormatByFileName(fname):
    import re
    if re.match(r"^.*\.fa(sta)?(\.gz)?$", fname):
        return "fasta"
    elif re.match(r"^.*\.((fq)|(fastq))(\.gz)?$", fname):
        return "fastq"
    else:
        raise "Unrecognized file type"

def error_profile(rcm, library, repertoire, limit=5):
    from Bio import SeqIO
    assert limit > 0
    with open(library) as f:
        reads = [rec for rec in SeqIO.parse(f, idFormatByFileName(library))]

    id2read = {str(rec.description): rec.seq for rec in reads}

    rcm = parse_rcm(rcm)

    with open(repertoire) as f:
        reads = [rec for rec in SeqIO.parse(f, "fasta")]

    cluster2center = {parse_cluster_mult(str(read.description))[0]: str(read.seq) for read in reads}

    from collections import defaultdict

    cluster2read = defaultdict(list)

    for id, cluster in rcm.iteritems():
        cluster = cluster
        if id in id2read:
            cluster2read[cluster].append(id2read[id])

    errors = defaultdict(int)
    errors01 = []
    nreads = 0
    errors_in_read = []
    for cluster, reads in cluster2read.iteritems():
        if len(reads) < limit:
            continue
        if cluster not in cluster2center:
            continue

        center = cluster2center[cluster]

        nreads += len(reads)
        for read in reads:
            er_in_read = 0
            for i in xrange(min(len(read), len(center)) - 21): # HACK!!!!
                if read[i] != center[i]:
                    errors[i] += 1
                    errors01.append(float(i) / (len(center) - 21))
                    er_in_read += 1
            errors_in_read.append(er_in_read)

    maxlen = max(errors.iterkeys())
    result = [0] * (maxlen + 1)
    for i, err in errors.iteritems():
        result[i] = err

    error_rate = float(len(errors01)) / nreads

    class Empty:
        pass

    res = Empty()
    res.error_profile = result
    res.errors01 = errors01
    res.error_rate = error_rate
    res.nreads = nreads
    res.limit = limit
    res.errors_in_read = errors_in_read
    return res








def rcm2Rand(filename="./igrc_out/final_repertoire.rcm"):
    import sklearn
    import sklearn.metrics
    import re

    reference = []
    clustering = []

    with open(filename) as fh:
        for line in fh:
            id, cluster = line.strip().split("\t")
            cluster = int(cluster)
            m = re.match("^antibody_(\\d+)_", id)
            ant = int(m.groups()[0])
            clustering.append(ant)
            reference.append(cluster)

    # find large >= clusteres
    from collections import defaultdict

    clust2mult = defaultdict(int)
    for c in clustering:
        clust2mult[c] += 1

    reference_large = []
    clustering_large = []
    for ref, cluster in zip(reference, clustering):
        if clust2mult[cluster] >= 5:
            reference_large.append(ref)
            clustering_large.append(cluster)

    # print reference
    # print clustering
    assert rand(clustering, clustering) == FM_index(clustering, clustering) == jaccard(clustering, clustering) == 1

    return FM_index(clustering, reference), rand(clustering, reference), rand_adj(clustering, reference), jaccard(clustering, reference), FM_index(clustering_large, reference_large), rand(clustering_large, reference_large), rand_adj(clustering_large, reference_large),jaccard(clustering_large, reference_large)


if __name__ == "__main__":
    print rcm2Rand()
