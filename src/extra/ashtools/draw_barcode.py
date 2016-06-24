#!/usr/bin/env python2

from Bio import SeqIO
import argparse
import sys
import Levenshtein
import igraph as ig
from align import align_reads
from align import consensus
from ig_basic import *
import numpy as np
import seaborn as sns
import os


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

def expand_path(path):
    from os.path import abspath, expanduser, expandvars

    return abspath(expanduser(expandvars(path)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot colored Hamming graph")
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument("--input", "-i",
                             nargs="+",
                             type=str,
                             help="Reads file(s) in FASTQ format, reads from one B-cell")
    input_group.add_argument("--file-list", "-I",
                             type=str,
                             help="Reads files list")
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
                        action="count",
                        help="Be verbose. Repeat several times for more verbosity")
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
                        default=1,
                        help="Minimal read length; shorter reads will be omited [default %(default)d]")
    parser.add_argument("-n", "--minimal-library-size",
                        type=int,
                        default=1,
                        help="Minimal library size; libs won't be processed [default %(default)d]")

    args = parser.parse_args()

    print("Reading subgraphs...")
    id2cluster = load_sgraphs(args.sgraphs)
    print("Subgraphs readed")

    if args.file_list is not None:
        with open(args.file_list, "r") as fh:
            args.input = [l.strip() for l in fh]

    print("The number of barcodes %d" % len(args.input))

    ftm_1 = np.zeros((4, 4), dtype=int)
    ftm_multinodes = np.zeros((4, 4), dtype=int)
    ftm_all = np.zeros((4, 4), dtype=int)

    for index in range(len(args.input)):
        fname = args.input[index]
        fname = expand_path(fname)

        templates = {"basename": os.path.basename(fname),
                     "absname": fname,
                     "dirname": os.path.dirname(fname),
                     "index": index}

        with open(fname, "rU") as fh:
            reads = list(SeqIO.parse(fh, "fastq"))

        reads = [read for read in reads if len(read.seq) >= args.minimal_read_length]


        lens = [len(read.seq) for read in reads]
        print "Min/max length: %d:%d" % (min(lens), max(lens))
        print sorted(lens)

        if len(reads) < args.minimal_library_size:
            print("Too few reads %d < %d" % (len(reads), args.minimal_library_size))
            continue
        print("%d reads using..." % len(reads))

        aligned_reads, max_shift = align_reads(reads, delta=args.delta, inplace=False)

        if args.delta != 0:
            print("Max shift %d" % max_shift)

        con = consensus(aligned_reads)
        if args.verbose > 2:
            print(con.ftm_all)
            print(con.ftm_1)
        ftm_1 += con.ftm_1
        ftm_all += con.ftm_all
        ftm_multinodes += con.ftm_multinodes

        lev_dists = [Levenshtein.distance(con.consensus, str(read.seq)) for read in aligned_reads]
        ham_dists = [Levenshtein.hamming(con.consensus, str(read.seq)) for read in aligned_reads]

        if args.drop_insertions:
            aligned_reads = [read for ld, hd, read in zip(lev_dists, ham_dists, aligned_reads) if hd == ld]
            print("Reads left %d" % len(aligned_reads))
            if len(aligned_reads) < 1:
                continue
            con = consensus(aligned_reads)

        str_reads = [str(read.seq) for read in aligned_reads]

        clusters = [id2cluster[read.id] if read.id in id2cluster else None for read in aligned_reads]

        main_cluster = most_popular_element(clusters)
        none_cluster_size = sum(cluster is None for cluster in clusters)
        main_cluster_size = sum(main_cluster == cluster for cluster in clusters)

        print("Main cluster %d/%d = %f" % (main_cluster_size,
                                        len(clusters),
                                        float(main_cluster_size) / len(clusters) if len(clusters) else float("nan")))
        print("Unclustered items %d" % none_cluster_size)
        print("Main cluster (excluded unclustered items) %d/%d = %f" % (main_cluster_size,
                                                                        len(clusters) - none_cluster_size,
                                                                        float(main_cluster_size) / (len(clusters) - none_cluster_size) if len(clusters) > none_cluster_size else float("nan")))

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
            pal = sns.color_palette("husl", n)
            pal = sns.color_palette("Set2")[:7]
            if len(pal) < n:
                # pal += sns.cubehelix_palette(n - len(pal), dark=0.5, light=0.95)
                pal += sns.color_palette("Blues_r", n - len(pal))
            colors_palette += [ig.drawing.colors.color_to_html_format(pal[i]) for i in range(n)]

        print("The number of different colors %d " % n_indices)
        print("The number of unique reads %d" % len(str_reads))

        if args.levenshtein:
            g = levenshtein_graph(str_reads, tau=args.tau)
        else:
            g = hamming_graph(str_reads, tau=args.tau, method="knuth")
            # g_naive = hamming_graph(str_reads, tau=args.tau, multiplicity=mult,
            #                         method="naive")
            # assert(sorted(g.get_edgelist()) == sorted(g_naive.get_edgelist()))

        if args.span and len(g.es):
            g = g.spanning_tree("weight")

        if len(g.es):
            g.es["label"] = map(int, g.es["weight"])
        g.vs["color"] = [colors_palette[i] for i in indices]
        g.vs["multiplicity"] = mult
        g.vs["label"] = map(int, g.vs["multiplicity"])
        # g.vs["label"] = indices

        vertex_sizes = [30 if m > 10 else 20 for m in mult]
        if len(g.es):
            edge_widths = [3 if w > 10 else 1 for w in g.es["weight"]]
            edge_colors = ["red" if w > 10 else "black" for w in g.es["weight"]]

        if args.figure:
            ig.plot(g, target=expand_path(args.figure % templates),
                    vertex_size=vertex_sizes,
                    edge_width=edge_widths,
                    edge_color=edge_colors,
                    bbox=(800, 800))

        if args.graph:
            g.write_dot(expand_path(args.graph % templates))

    print("Seq errors")
    print(ftm_1)
    print("All errors")
    print(ftm_all)
    print("Amplification (multinodes) errors")
    print(ftm_multinodes)
