#!/usr/bin/env python2

from Bio import SeqIO
import argparse
import hashlib


from align import align_reads


def get_yanas_size(s):
    import re
    m = re.match(r".*___(\d+)$", s)
    return int(m.groups()[0])

assert(get_yanas_size("cluster___217733___size___77") == 77)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Construct HG for repertoire")
    parser.add_argument("input",
                        type=str,
                        help="FASTA file with repertoire")
    parser.add_argument("output",
                        type=str,
                        help="output prefix")
    parser.add_argument("--rep-out",
                        type=str,
                        help="FASTA file for resulting repertoire")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Be verbose (default no)")
    # parser.add_argument("--igrc",
    #                     action="store_true",
    #                     help="Whether repertoire obtained from IG_RC")
    parser.add_argument("--make-graphs", "-g",
                        action="store_true",
                        help="Whether to make graphs")
    parser.add_argument("-d", "--delta",
                        type=int,
                        default=0,
                        help="Max possible shift from reference (the largest string) [default %(default)d]")

    args = parser.parse_args()

    with open(args.input, "rU") as fh:
        reads = list(SeqIO.parse(fh, "fasta"))

    reads = [read for read in reads if len(read.seq) >= 100]


    # if args.igrc:
    reads = [read for read in reads if get_yanas_size(read.id) >= 50]

    print("%d reads using..." % len(reads))

    aligned_reads, max_shift = align_reads(reads, delta=args.delta, inplace=False)

    from ig_basic import *

    def save_span(g, fname=None):
        import re
        import igraph as ig
        tree = g.spanning_tree("weight")
        tree.es["label"] = map(int, tree.es["weight"])
        tree.vs["label"] = map(int, tree.vs["multiplicity"])
        ig.plot(tree, target=fname, vertex_size=25, bbox=(1800, 1800))
        tree.write_dot(re.sub(r"\.[^.]+$", ".dot", fname))

    str_reads = [str(read.seq) for read in aligned_reads]
    str_reads, mult = count_multiplicity(str_reads)

    print("The number of unique reads %d" % len(str_reads))

    if args.rep_out is None:
        args.rep_out = args.output + ".fasta"

    if args.rep_out is not None:
        with open(args.rep_out, "w") as fh:
            for read in aligned_reads:
                fh.write(">%s\n%s\n" % (hashlib.sha224(str(read.seq)).hexdigest(),
                                        read.seq))

    if args.make_graphs:
        g = hamming_graph(str_reads, tau=400000, multiplicity=mult)
        save_span(g, "%s_rep_hamming_stripped.pdf" % args.output)

        g = levenshtein_graph(str_reads, tau=400000, multiplicity=mult)
        save_span(g, "%s_rep_leven_stripped.pdf" % args.output)

    #
    # str_reads = [str(read.seq) for read in reads]
    # str_reads, mult = count_multiplicity(str_reads)
    # g = levenshtein_graph(str_reads, tau=400000, multiplicity=mult)
    # save_span(g, "%s_rep_leven.pdf" % args.output)
