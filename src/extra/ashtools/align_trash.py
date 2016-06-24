sys.exit(0)

import numpy as np
import Levenshtein

lev_dists = [Levenshtein.distance(con.consensus, str(read.seq)) for read in aligned_reads]
ham_dists = [Levenshtein.hamming(con.consensus, str(read.seq)) for read in aligned_reads]


filtered_reads = [read for ld, hd, read in zip(lev_dists, ham_dists, aligned_reads) if hd == ld]
if not args.filtered is None:
    with smart_open(args.filtered) as fh:
        SeqIO.write(filtered_reads, fh, "fastq")

con = consensus(filtered_reads)

print("Max Lev %d, max Ham %d" % (max(lev_dists), max(ham_dists)))
print(zip(lev_dists, ham_dists))

most_mismatched_read = aligned_reads[np.argmax(con.reads_mismatches)]
print(con.consensus)
print(most_mismatched_read.seq)
print(con.mismatches)
# print(con.first_votes)
# print(con.second_votes)
# print(con.first_votes - con.second_votes)
print(con.reads_mismatches)
print("Max mismatches in one read %d" % max(con.reads_mismatches))


from ig_basic import *
import igraph as ig
def save_span(g, fname=None):
    tree = g.spanning_tree("weight")
    tree.es["label"] = map(int, tree.es["weight"])
    tree.vs["label"] = map(int, tree.vs["multiplicity"])
    ig.plot(tree, target=fname, vertex_size=25, bbox=(1800, 1800))

# str_reads = [str(read.seq) for read in reads]
str_reads = [str(read.seq) for read in filtered_reads]
str_reads, mult = count_multiplicity(str_reads)
g = hamming_graph(str_reads, tau=400000, multiplicity=mult)
save_span(g, "hamming_filtered.pdf")


str_reads = [str(read.seq) for read in aligned_reads]
str_reads, mult = count_multiplicity(str_reads)
g = hamming_graph(str_reads, tau=400000, multiplicity=mult)
save_span(g, "hamming_stripped.pdf")


str_reads = [str(read.seq) for read in aligned_reads]
str_reads, mult = count_multiplicity(str_reads)
g = levenshtein_graph(str_reads, tau=400000, multiplicity=mult)
save_span(g, "leven_stripped.pdf")

