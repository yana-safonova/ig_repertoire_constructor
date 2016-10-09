#!/usr/bin/env python2

from argparse import ArgumentParser
from Bio import SeqIO
import sys
import os
import os.path
from collections import defaultdict

import numpy as np

def get_data(filename):
    CDR3s = []
    mutationCounts = []

    with open(filename) as f:
        header = next(f)
        for line in f:
            items = line.strip().split("\t")
            CDR3 = items[2]
            CDR3s.append(CDR3)
            bestVAlignment = items[5]
            targetFrom, targetTo, targetLength, queryFrom, queryTo, mutations, alignmentScore = bestVAlignment.strip().split("|")
            substitutions = mutations.count("S")
            insertions = mutations.count("I")
            deletions = mutations.count("D")

            mutationCount = substitutions + insertions + deletions
            mutationCounts.append(mutationCount)

    return CDR3s, mutationCounts


class DSU:
    def __init__(self, n):
        self.data = range(n)

    def get(self, i):
        if self.data[i] == i:
            return i
        else:
            self.data[i] = self.get(self.data[i])
            return self.data[i]

    def join(self, i, j):
        i = self.get(i)
        j = self.get(j)

        assert self.data[i] == i
        assert self.data[j] == j
        self.data[i] = j

    def getAll(self):
        for i in xrange(len(self.data)):
            self.get(i)
        return self.data


def hamming(s1, s2):
    if len(s1) != len(s2):
        return 10005000

    res = 0
    for l1, l2 in zip(s1, s2):
        if l1 != l2:
            res += 1

    return res

import itertools


def write_cdr3s(cdr3s, filename):
    with open(filename, "w") as f:
        for i, s in enumerate(cdr3s):
            f.write(">id%d\n%s\n" % (i, s))


current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = current_dir + "/../../"


def call_swg(inputfile, outputfile):
    run_graph_constructor = igrec_dir + "/build/release/bin/ig_swgraph_construct"
    cmd = run_graph_constructor + " -i " + inputfile + " -o " + outputfile + " --tau 3 -T 0 -S1 -k8"
    os.system(cmd)


def parse_conn_comp_ids(graph_file):
    with open(graph_file) as g:
        header = next(g)
        V, E, _ = header.strip().split()
        V = int(V)
        E = int(E)

        dsu = DSU(V)

        for i, line in enumerate(g):
            line = line.strip().split()
            neibs = line[0::2]
            neibs = map(int, neibs)
            neibs = [j - 1 for j in neibs]
            for j in neibs:
                dsu.join(i, j)

    return dsu.getAll()


def conn_comp(cdrs):
    import tempfile

    CDR_fa = tempfile.mkstemp(suffix=".fa", prefix="nanodivan_")[1]
    CDR_graph = tempfile.mkstemp(suffix=".graph", prefix="nanodivan_")[1]
    write_cdr3s(cdrs, CDR_fa)
    call_swg(CDR_fa, CDR_graph)
    ids = parse_conn_comp_ids(CDR_graph)

    os.remove(CDR_fa)
    os.remove(CDR_graph)

    return ids


if __name__ == "__main__":
    filename = "/ssd/ig_repertoire_constructor/py/nano_divan/test7_mixcr2/features.txt"
    filename = sys.argv[1]

    CDR3s, mutationCounts = get_data(filename)
    # print CDR3s
    print "Mean mutaion count in V segment:", np.mean(mutationCounts)

    abundances = defaultdict(int)
    for cdr3 in CDR3s:
        abundances[cdr3] += 1

    unique_cdr3s = abundances.keys()
    abundances = [abundances[key] for key in unique_cdr3s]
    print "# unique CDR3", len(abundances)
    print "Max CDR3 abundance", max(abundances)

    Ns = np.array(abundances, dtype=float)
    freqs = Ns / sum(Ns)

    shannon_index = -sum(freqs * np.log(freqs))
    simpson_index = sum(freqs ** 2)

    print "Shannon:", shannon_index
    print "Simpson:", simpson_index

    # dsu = DSU(len(unique_cdr3s))
    #
    #
    # write_cdr3s(unique_cdr3s, "CDR.fa")
    # call_swg("CDR.fa", "CDR.graph")
    #
    # tau = 3
    # for i, j in itertools.combinations(range(len(unique_cdr3s)), 2):
    #     if hamming(unique_cdr3s[i], unique_cdr3s[j]) <= tau:
    #         dsu.join(i, j)
    # labels = dsu.getAll()
    #
    labels = conn_comp(unique_cdr3s)

    compsizecdrs = defaultdict(int)
    compsizeuniquecdrs = defaultdict(int)

    for label, size in zip(labels, Ns):  # TODO check consistensy keys-values
        compsizecdrs[label] += size
        compsizeuniquecdrs[label] += 1

    Nsignletons = sum(1 for size in compsizeuniquecdrs.values() if size == 1)
    print "Singletons: ", Nsignletons

    CNs = np.array(compsizecdrs.values(), dtype=float)
    assert sum(CNs) == sum(Ns)
    freqs = CNs / sum(CNs)

    shannon_index = -sum(freqs * np.log(freqs))
    simpson_index = sum(freqs ** 2)

    print "Shannon:", shannon_index
    print "Simpson:", simpson_index
