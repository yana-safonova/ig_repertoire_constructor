def lfactorial(x):
    """
    Compute log-factorial using log-gamma function
    """
    from math import lgamma

    return lgamma(x + 1)


def anonnamedtuple(field_names, *args, **kwargs):
    """
    Make anonymous named tuple (like make_tuple in C++11)
    """
    from collections import namedtuple

    return namedtuple("_", field_names)(*args, **kwargs)


def groupby_dict(iterable, key=None, lazy=False, key_as_index=False,
                 raw_dict=False):
    """
    Analogue of itertools.groupby, but using dict instead of sorted list
    """
    from collections import defaultdict

    if key is None:
        key = lambda x: x

    d = defaultdict(list)
    if not key_as_index:
        for e in iterable:
            d[key(e)].append(e)
    else:
        for k, e in zip(key, iterable):
            d[k].append(e)

    if raw_dict:
        return d
    else:
        return d.iteritems() if lazy else list(d.iteritems())


def count_multiplicity(reads):
    gd = groupby_dict(reads, reads, key_as_index=True)
    reads = []
    multiplicity = []
    for read, l in gd:
        reads.append(read)
        multiplicity.append(len(l))

    return reads, multiplicity


def hamming_graph_knuth(reads, tau=1, **kwargs):
    """
    Construct hamming(tau) graph using Knuth's algorithm
    """
    import igraph as ig
    from collections import defaultdict
    from Levenshtein import hamming

    l = len(reads[0])

    for read in reads:
        assert(len(read) == l)

    piece_len = int(l / (tau + 1))
    piece_len_last = l - piece_len * tau

    piece_lens = [piece_len] * tau + [piece_len_last]

    g = ig.Graph(len(reads))
    g.vs["read"] = reads

    for attr_name, attr_data in kwargs.iteritems():
        g.vs[attr_name] = attr_data

    edges_for_check = set()

    for j in range(tau + 1):
        sets = defaultdict(list)
        for i in range(len(reads)):
            substr = reads[i][j * piece_len:j * piece_len + piece_lens[j]]
            sets[substr].append(i)

        for v_list in sets.itervalues():
            # print("N += %d" % len(v_list))
            for i1 in range(len(v_list)):
                for i2 in range(i1 + 1, len(v_list)):
                    edges_for_check.add((v_list[i1], v_list[i2]))

    # print("Edges for check %d" % len(edges_for_check))
    for v1, v2 in edges_for_check:
        read1, read2 = reads[v1], reads[v2]
        d = hamming(read1, read2)
        if d <= tau:
            g.add_edge(v1, v2, weight=d)

    return g


def hamming_graph_vp(reads, tau=1, **kwargs):
    """
    Construct hamming(tau) graph using Vantage point tree
    """
    import igraph as ig
    from Levenshtein import hamming
    from vptree import VPTree

    l = len(reads[0])

    for read in reads:
        assert(len(read) == l)

    g = ig.Graph(len(reads))
    g.vs["read"] = reads

    for attr_name, attr_data in kwargs.iteritems():
        g.vs[attr_name] = attr_data

    def dist(i, j):
        return hamming(reads[i], reads[j])

    tree = VPTree(range(len(reads)), dist)

    for i in range(len(reads)):
        neib = tree.get_all_in_range(i, tau)
        for d, j in neib:
            assert(d <= tau)
            if j > i:
                g.add_edge(i, j, weight=d)

    return g


def hamming_graph_bk(reads, tau=1, **kwargs):
    """
    Construct hamming(tau) graph using Vantage point tree
    """
    import igraph as ig
    from Levenshtein import hamming
    from vptree import BKTree

    l = len(reads[0])

    for read in reads:
        assert(len(read) == l)

    g = ig.Graph(len(reads))
    g.vs["read"] = reads

    for attr_name, attr_data in kwargs.iteritems():
        g.vs[attr_name] = attr_data

    def dist(i, j):
        return hamming(reads[i], reads[j])

    tree = BKTree(range(len(reads)), dist)

    for i in range(len(reads)):
        neib = tree.get_all_in_range(i, tau)
        for d, j in neib:
            assert(d <= tau)
            if j > i:
                g.add_edge(i, j, weight=d)

    return g


def hamming_graph_naive(reads, tau=1, **kwargs):
    """
    Construct hamming(tau) graph using naive O(N**2 d) algorithm
    """
    import igraph as ig
    import numpy as np
    from Levenshtein import hamming

    l = len(reads[0])

    for read in reads:
        assert(len(read) == l)

    N = len(reads)
    m = np.zeros((N, N), dtype=int)

    for i in range(N):
        for j in range(i):
            dist = hamming(reads[i], reads[j])
            m[i, j] = m[j, i] = dist if dist <= tau else 0

    # Be careful! Zero elements are not interpreted as zero-length edges
    g = ig.Graph.Weighted_Adjacency(m.tolist(),
                                    mode="UNDIRECTED",
                                    attr="weight",
                                    loops=False)

    g.vs["read"] = reads

    for attr_name, attr_data in kwargs.iteritems():
        g.vs[attr_name] = attr_data

    return g


def hamming_graph(reads, tau=1, method="knuth", **kwargs):
    """
    Hamming graph(tau) construction wrapper
    """
    l = len(reads[0])

    for read in reads:
        assert(len(read) == l)

    if l / (tau + 1) < 3 and method == "knuth":
        method = "naive"

    if method == "naive":
        g = hamming_graph_naive(reads, tau=tau, **kwargs)
    elif method == "knuth":
        g = hamming_graph_knuth(reads, tau=tau, **kwargs)
    elif method == "vp":
        g = hamming_graph_vp(reads, tau=tau, **kwargs)
    elif method == "bk":
        g = hamming_graph_bk(reads, tau=tau, **kwargs)
    else:
        raise ValueError("unknown `method`")

    return g


def levenshtein_graph(reads, tau=1, **kwargs):
    """
    Construct Levenshtein(tau) graph using naive O(N**2 d) algorithm
    """
    import igraph as ig
    import numpy as np
    from Levenshtein import distance as levenshtein_dist

    N = len(reads)
    m = np.zeros((N, N), dtype=int)

    for i in range(N):
        for j in range(i):
            dist = levenshtein_dist(reads[i], reads[j])
            m[i, j] = m[j, i] = dist if dist <= tau else 0

    # Be careful! Zero elements are not interpreted as zero-length edges
    g = ig.Graph.Weighted_Adjacency(m.tolist(),
                                    mode="UNDIRECTED",
                                    attr="weight",
                                    loops=False)

    g.vs["read"] = reads

    for attr_name, attr_data in kwargs.iteritems():
        g.vs[attr_name] = attr_data

    return g


def findDuplicates(l):
    """
    Find and count duplicates in list using dict
    """
    from collections import defaultdict, namedtuple

    element_count = defaultdict(int)
    for e in l:
        element_count[e] += 1

    uniques, counts = [list(_) for _ in zip(*element_count.iteritems())]

    return namedtuple("elements_counts",
                      ["elements", "counts"])(uniques, counts)


def count_triplet(s, t, align=0):
    return sum([s[i:(i + 3)] == t for i in range(align, len(s) - 3 + 1, 3)])


def stop_codons_count(s, align=0):
    """
    Stop codons (from wiki):
        - TAG ("amber")
        - TAA ("ochre")
        - TGA ("opal" or "umber")
    """
    cs = [count_triplet(s, t, align=align) for t in ["TAG", "TAA", "TGA"]]
    return sum(cs)


def optimal_align_using_stop_codons(reads):
    import numpy as np
    return np.argmin([sum([stop_codons_count(read, align=align)
                           for read in reads])
                      for align in [0, 1, 2]])


def stop_codon_indices(s):
    res = []
    for i in range(len(s) - 3 + 1):
        if s[i:i + 3] in ["TAG", "TAA", "TGA"]:
            res.append(i)
    return res


def stop_codon_indices_all(reads):
    res = []
    for read in reads:
        res.extend(stop_codon_indices(read))

    return res


def extract_barcode_age(s):
    """
    @M01691:10:000000000-A7F7L:1:2107:21387:4106_2:N:0:1_UMI:CTGACGTTACTCGG:GGGGGGGGGGGGGG_517009
    ->
    CTGACGTTACTCGG
    """
    import re

    m = re.match(r".*UMI:([ACTG]{12,17}):.*", s)
    return m.groups()[0] if m else None


def extract_barcode_yale(s):
    """
    1_merged_read_MISEQ:166:000000000-AAYHL:1:1101:11497:1414_1:N:0:TGACCA|SEQORIENT=F|PRIMER=Human-IGK|BARCODE=ACGTGTAGTAAAATAAT
    ->
    A CGTG TAGT AAAA TAAT
    """
    import re

    m = re.match(r".*BARCODE=([ACTG]{17}).*", s)
    return m.groups()[0] if m else None


def extract_barcode(s):
    bc = extract_barcode_age(s)
    if bc is not None:
        return bc
    else:
        return extract_barcode_yale(s)

assert extract_barcode("2_merged_read_MISEQ:166:000000000-AAYHL:1:1101:19266:1412_1:N:0:TGACCA|SEQORIENT=F|PRIMER=Human-IGK|BARCODE=CTGCGAGCACGACTACG") == "CTGCGAGCACGACTACG"
assert extract_barcode("1_merged_read_MISEQ:166:000000000-AAYHL:1:1101:11497:1414_1:N:0:TGACCA|SEQORIENT=F|PRIMER=Human-IGK|BARCODE=ACGTGTAGTAAAATAAT") == "ACGTGTAGTAAAATAAT"

assert extract_barcode("@M01691:10:000000000-A7F7L:1:2107:21387:4106_2:N:0:1_UMI:CTGACGTTACTCGG:GGGGGGGGGGGGGG_517009") == "CTGACGTTACTCGG"
assert extract_barcode("@M01691:10:000000000-A7F7L:1:2109:5441:12945_1:N:0:1_UMI:AAGTATTTAGTAAC:FFEFFGGGGGGFG9_565461") == "AAGTATTTAGTAAC"
assert extract_barcode("M01691:56:000000000-ABFHP:1:1101:11306:1728 1:N:0:1 UMI:CTTTACGTACGT:<FFGGGGFF@@F") == "CTTTACGTACGT"

assert extract_barcode("@98_merged_read_HWI-M01418:310:000000000-AAH58:1:1101:22293:2129_1:N:0:ATTCCT_UMI:TATAGAATATAATTTAA:-AC<,,C,,,<CF<,CD") == "TATAGAATATAATTTAA"


def most_popular_element(l):
    from collections import defaultdict

    count = defaultdict(int)
    for e in l:
        if e is not None:
            count[e] += 1

    if len(count) == 0:
        return None
    else:
        return max(count.keys(), key=lambda k: count[k])
