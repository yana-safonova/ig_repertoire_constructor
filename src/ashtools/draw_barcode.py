#!/usr/bin/env python2

from Bio import SeqIO
import argparse
import sys
import Levenshtein
import igraph as ig
from align import align_reads
from align import consensus
from ig_basic import *


def load_sgraphs(dir):
    import glob
    import os

    dir = os.path.abspath(dir)

    id2cluster = {}

    for decomposition in glob.glob("%s/dense_subgraphs/*.txt" % dir):
        basename = os.path.basename(decomposition)
        with open(decomposition, "r") as dec_fh:
            for line in dec_fh:
                cluster, seq_id, align = line.split()
                cluster = int(cluster)
                seq_id = seq_id.strip()
                id2cluster[seq_id] = (basename, cluster)

    return id2cluster


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
    parser.add_argument("-U", "--omit-unclustered",
                        action="store_true",
                        help="Whether to omit reads which was unclustered by igRC")
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

    if args.delta != 0:
        print("Max shift %d" % max_shift)

    con = consensus(aligned_reads)

    lev_dists = [Levenshtein.distance(con.consensus, str(read.seq)) for read in aligned_reads]
    ham_dists = [Levenshtein.hamming(con.consensus, str(read.seq)) for read in aligned_reads]

    if args.drop_insertions:
        aligned_reads = [read for ld, hd, read in zip(lev_dists, ham_dists, aligned_reads) if hd == ld]
        con = consensus(aligned_reads)
        print("Reads left %d" % len(aligned_reads))

    str_reads = [str(read.seq) for read in aligned_reads]
    id2cluster = load_sgraphs(args.sgraphs)
    clusters = [id2cluster[read.id] if read.id in id2cluster else None for read in aligned_reads]

    if args.omit_unclustered:
        _ = [(str_read, cluster) for str_read, cluster in zip(str_reads, clusters) if cluster is not None]
        str_reads, clusters = map(list, zip(*_))

    gd = groupby_dict(clusters, str_reads,
                      lazy=False, key_as_index=True, raw_dict=True)

    str_reads, mult = count_multiplicity(str_reads)
    clusters = [most_popular_element(gd[str_read]) for str_read in str_reads]

    clusters_not_none = [cluster for cluster in clusters if cluster is not None]
    clusters_unique, clusters_mult = count_multiplicity(clusters_not_none)
    clusters_unique = [_[1] for _ in sorted(zip(clusters_mult,
                                                clusters_unique),
                                            reverse=True)]

    read2index = {cluster: i for cluster, i in zip(clusters_unique,
                                                   range(1, len(clusters_unique) + 1))}

    read2index[None] = 0

    indices = [read2index[cluster] for cluster in clusters]

    n_indices = len(set(indices))
    colors_palette = ["white", "red"]

    if len(colors_palette) < max(indices) + 1:
        n = max(indices) + 1 - len(colors_palette)
        pal = ig.drawing.colors.RainbowPalette(n=n, start=0.3)
        colors_palette += [ig.drawing.colors.color_to_html_format(pal.get(i)) for i in range(n)]

    print("The number of different colors %d " % n_indices)
    print("The number of unique reads %d" % len(str_reads))

    if args.levenshtein:
        g = levenshtein_graph(str_reads, tau=args.tau, multiplicity=mult)
    else:
        g = hamming_graph(str_reads, tau=args.tau, multiplicity=mult)

    if args.span:
        g = g.spanning_tree("weight")

    g.es["label"] = map(int, g.es["weight"])
    g.vs["color"] = [colors_palette[i] for i in indices]
    g.vs["label"] = map(int, g.vs["multiplicity"])

    if args.figure:
        ig.plot(g, target=args.figure,
                vertex_size=25, bbox=(1800, 1800))

    if args.graph:
        g.write_dot(args.graph)
