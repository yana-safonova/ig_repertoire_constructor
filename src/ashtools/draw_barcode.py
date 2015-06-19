#!/usr/bin/env python2

from Bio import SeqIO
import argparse
import sys
import Levenshtein
import igraph as ig


def unique_counter_filler(l):
    m = {}
    res = []
    i = 0
    for e in l:
        if e not in m:
            i = i + 1
            m[e] = ci = i
        else:
            ci = m[e]
        res.append(ci)

    return res


def load_sgraphs(dir):
    import glob
    import os

    dir = os.path.abspath(dir)

    id2cluster = {}

    dec_index = 0
    for decomposition in glob.glob("%s/dense_subgraphs/*.txt" % dir):
        dec_index += 1
        with open(decomposition, "r") as dec_fh:
            for line in dec_fh:
                cluster, seq_id, align = line.split()
                cluster = int(cluster)
                seq_id = seq_id.strip()
                id2cluster[seq_id] = (dec_index, cluster)

    return id2cluster

# from read_barcode_splitter import extract_barcode
from align import align_reads


def get_colors(ids, id2cluster):
    import igraph as ig

    idx = [id2cluster[_] if _ in id2cluster else 0 for _ in ids]
    idx = unique_counter_filler(idx)

    colors = ["black", "red", "blue", "pink", "violet", "grey", "magenta", "yellow", "green", "green4", "green", "green"] + ig.drawing.colors.known_colors.keys()

    return [colors[i] for i in idx]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot colored Hamming graph")
    parser.add_argument("input",
                        type=str,
                        help="Reads file in FASTQ format, reads from one B-cell")
    parser.add_argument("--sgraphs",
                        type=str,
                        help="Dir with sgraphs")
    parser.add_argument("--figure", "-F",
                        type=str,
                        help="File for graph image output")
    parser.add_argument("--graph", "-G",
                        type=str,
                        help="File for graph output in DOT format")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Be verbose (default no)")
    parser.add_argument("-t", "--tau",
                        type=int,
                        default=3,
                        help="Hamming graph maximal length [default %(default)d]")
    parser.add_argument("-S", "--span",
                        action="store_true",
                        help="Whether to extract spanning tree before graph plotting")
    parser.add_argument("-L", "--levenshtein",
                        action="store_true",
                        help="Whether to use Levenshtein distance instead of Hamming")
    parser.add_argument("-D", "--drop-insertions",
                        action="store_true",
                        help="Whether to drop reads with insertions")
    parser.add_argument("-d", "--delta",
                        type=int,
                        default=0,
                        help="Max possible shift from reference (the largest string) [default %(default)d]")
    parser.add_argument("-l", "--minimal-read-length",
                        type=int,
                        default=0,
                        help="Minimal read length; shorter reads will be omited [default %(default)d]")
    parser.add_argument("-n", "--minimal-library-size",
                        type=int,
                        default=0,
                        help="Minimal library size; libs won't be processed [default %(default)d]")

    args = parser.parse_args()
    with open(args.input, "rU") as fh:
        reads = list(SeqIO.parse(fh, "fastq"))

    reads = [read for read in reads if len(read.seq) >= args.minimal_read_length]

    if len(reads) < args.minimal_library_size:
        print("Too few reads %d < %d" % (len(reads), args.minimal_library_size))
        sys.exit(0)
    print("%d reads using..." % len(reads))

    aligned_reads, max_shift = align_reads(reads, delta=args.delta, inplace=False)

    from align import consensus
    from ig_basic import count_multiplicity

    con = consensus(aligned_reads)

    lev_dists = [Levenshtein.distance(con.consensus, str(read.seq)) for read in aligned_reads]
    ham_dists = [Levenshtein.hamming(con.consensus, str(read.seq)) for read in aligned_reads]

    if args.drop_insertions:
        aligned_reads = [read for ld, hd, read in zip(lev_dists, ham_dists, aligned_reads) if hd == ld]
        con = consensus(aligned_reads)
        print("Reads left %d" % len(aligned_reads))

    str_reads = [str(read.seq) for read in aligned_reads]
    id2cluster = load_sgraphs(args.sgraphs)
    # print(id2cluster.keys())

    colors = get_colors([read.id for read in aligned_reads], id2cluster)


    read2color = {}
    for read, color in zip(aligned_reads, colors):
        read2color[str(read.seq)] = color

    str_reads, mult = count_multiplicity(str_reads)

    print("The number of unique reads %d" % len(str_reads))

    from ig_basic import *

    if args.levenshtein:
        g = levenshtein_graph(str_reads, tau=args.tau, multiplicity=mult)
    else:
        g = hamming_graph(str_reads, tau=args.tau, multiplicity=mult)

    if args.span:
        g = g.spanning_tree("weight")

    g.es["label"] = map(int, g.es["weight"])
    g.vs["color"] = [read2color[read] for read in g.vs["read"]]
    g.vs["label"] = map(int, g.vs["multiplicity"])

    if args.figure:
        ig.plot(g, target=args.figure,
                vertex_size=25, bbox=(1800, 1800))

    if args.graph:
        g.write_dot(args.graph)
