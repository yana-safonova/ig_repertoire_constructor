#!/usr/bin/env python2

from Bio import SeqIO
import argparse
import sys
import Levenshtein
import numpy as np


def count_matches(s, reference, align=0):
    if align < 0:
        s, reference = reference, s
        align = -align

    reference, s = str(reference), str(s)

    return sum(c1 == c2 for c1, c2 in zip(s, reference[align:]))


def optimal_align(s, reference, delta=10):
    return max(range(-delta, delta + 1),
               key=lambda align: count_matches(s, reference, align))


def align_reads(reads, delta=19, inplace=False):
    # Use the longest read as "reference"
    reference = max(reads, key=lambda read: len(read.seq))

    aligns = [optimal_align(read.seq, reference.seq, delta=delta) for read in reads]

    max_align = max(aligns)
    aligns = [align - max_align for align in aligns]

    if not inplace:
        reads = reads[:]

    reads = [read[-align:] for read, align in zip(reads, aligns)]

    min_len = min([len(read.seq) for read in reads])

    reads = [read[:min_len] for read in reads]

    return reads, -min(aligns)


from collections import namedtuple

Consensus = namedtuple("Consensus",
                       "consensus mismatches reads_mismatches first_votes second_votes ftm_all ftm_1 ftm_multinodes")


def consensus(reads, raw_strings=False):
    import numpy as np

    if not raw_strings:
        reads = [str(read.seq) for read in reads]

    n = max(map(len, reads))
    for read in reads:
        assert(len(read) == n)

    nucleotides = np.array(['A', 'C', 'G', 'T'])
    nuc_ind = {nucleotides[i]: i for i in range(len(nucleotides))}

    consm = np.zeros((4, n), dtype=int)
    for j in range(len(reads)):
        s = reads[j]
        for i in range(n):
            consm[nuc_ind[s[i]], i] += 1

    consensus = "".join(nucleotides[np.argmax(consm, axis=0)])

    mismatches = np.sum(consm, axis=0) - np.max(consm, axis=0)

    first_votes = np.max(consm, axis=0)
    second_votes = np.percentile(consm, q=66, axis=0, interpolation="nearest")
    trird_votes = np.percentile(consm, q=33, axis=0, interpolation="nearest")
    fourth_votes = np.min(consm, axis=0)

    reads_mismatches = [n - count_matches(read, consensus, align=0) for read in reads]

    ftm_all = np.zeros((4, 4), dtype=int)
    ftm_1 = np.zeros((4, 4), dtype=int)
    ftm_multinodes = np.zeros((4, 4), dtype=int)
    multinode_threshold = 2

    for read in reads:
        for _from, _to in zip(consensus, read):
            ftm_all[nuc_ind[_from], nuc_ind[_to]] += 1

    for read, mm in zip(reads, reads_mismatches):
        if mm <= 1:
            for _from, _to in zip(consensus, read):
                ftm_1[nuc_ind[_from], nuc_ind[_to]] += 1

    from ig_basic import count_multiplicity
    unique_reads, mult = count_multiplicity(reads)
    for read, mul in zip(unique_reads, mult):
        if mul >= multinode_threshold:
            for _from, _to in zip(consensus, read):
                ftm_multinodes[nuc_ind[_from], nuc_ind[_to]] += mul

    return Consensus(consensus=consensus,
                     mismatches=mismatches,
                     first_votes=first_votes,
                     second_votes=second_votes,
                     reads_mismatches=reads_mismatches,
                     ftm_all=ftm_all,
                     ftm_1=ftm_1,
                     ftm_multinodes=ftm_multinodes)


import contextlib
@contextlib.contextmanager
def smart_open(filename=None):
    """
    From http://stackoverflow.com/questions/17602878/how-to-handle-both-with-open-and-sys-stdout-nicely
    """
    import sys

    if filename and filename != '-':
        fh = open(filename, 'w')
    else:
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Align reads to each other usign naive approach")
    parser.add_argument("input",
                        type=str,
                        help="Reads file in FASTQ format")
    parser.add_argument("output",
                        type=str,
                        help="Output FASTQ file, \"-\" for stdout")
    parser.add_argument("--filtered",
                        type=str,
                        help="Output FASTQ file for filtered reads, \"-\" for stdout")
    parser.add_argument("--consensus",
                        type=str,
                        default="consensus.fasta",
                        help="Output FASTA file consensus")
    parser.add_argument("--stats",
                        type=str,
                        help="Output file for stats")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Be verbose (default no)")
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

    n_initial = len(reads)

    reads = [read for read in reads if len(read.seq) >= args.minimal_read_length]

    if len(reads) < args.minimal_library_size:
        print("Too few reads %d < %d" % (len(reads), args.minimal_library_size))
        sys.exit(0)

    aligned_reads, max_shift = align_reads(reads, delta=args.delta, inplace=False)

    if args.verbose:
        print("Maximal shift: %d" % max_shift)
        print("Seq length: %d" % len(aligned_reads[0].seq))
        print("The number of reads: %d" % len(aligned_reads))

    with smart_open(args.output) as fh:
        SeqIO.write(aligned_reads, fh, "fastq")

    con = consensus(aligned_reads)

    if args.consensus is not None:
        with open(args.consensus, "a") as fh:
            fh.write(">%s___size___%d\n%s\n" % (args.input,
                                                len(aligned_reads),
                                                con.consensus))


    lev_dists = [Levenshtein.distance(con.consensus, str(read.seq)) for read in aligned_reads]
    ham_dists = [Levenshtein.hamming(con.consensus, str(read.seq)) for read in aligned_reads]

    filtered_reads = [read for ld, hd, read in zip(lev_dists, ham_dists, aligned_reads) if hd == ld]
    if not args.filtered is None:
        with smart_open(args.filtered) as fh:
            SeqIO.write(filtered_reads, fh, "fastq")

    con_filtered = consensus(filtered_reads)
    ham_dists_filtered = [Levenshtein.hamming(con_filtered.consensus,
                                              str(read.seq)) for read in filtered_reads]

    # Compute some stats
    n_initial = n_initial
    n_all = len(reads)
    n_filtered = len(filtered_reads)
    max_mismatched_read_dist = max(con.reads_mismatches)

    ham_dists = np.array(ham_dists)
    ham_dists_filtered = np.array(ham_dists_filtered)
    lev_dists = np.array(lev_dists)

    neigh_count = []
    for d in [0, 1, 2, 3, 4, 5, 6, 7]:
        neigh_count.append(np.sum(ham_dists == d))

    neigh_count_f = []
    for d in [0, 1, 2, 3, 4, 5, 6, 7]:
        neigh_count_f.append(np.sum(ham_dists_filtered == d))

    pois_lam_filtered = np.mean(ham_dists_filtered)

    template = "%s" + " %d" * 3 + " " + " %d" * 8 + " " + " %d" * 8 + "  %f\n"
    if args.stats is not None:
        with open(args.stats, "a") as fh:
            fh.write(template %
                     (args.input,
                      n_initial, n_all, n_filtered,
                      neigh_count[0],
                      neigh_count[1],
                      neigh_count[2],
                      neigh_count[3],
                      neigh_count[4],
                      neigh_count[5],
                      neigh_count[6],
                      neigh_count[7],
                      neigh_count_f[0],
                      neigh_count_f[1],
                      neigh_count_f[2],
                      neigh_count_f[3],
                      neigh_count_f[4],
                      neigh_count_f[5],
                      neigh_count_f[6],
                      neigh_count_f[7],
                      pois_lam_filtered))
