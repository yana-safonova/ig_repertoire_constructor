#!/usr/bin/env python2

from argparse import ArgumentParser
from Bio import SeqIO
import sys
import os.path
import tempfile
from collections import defaultdict
from ig_remove_low_abundance_reads import smart_open
from rcm_utils import read_rcm_list, read_rcm_map, rcm2rcmint, combine_rcms


def parse_cluster_mult(id):
    import re
    id = str(id)
    m = re.match(r"^cluster___(\d+)___size___(\d+)$", id)
    if m:
        g = m.groups()
        cluster = int(g[0])
        mult = int(g[1])
        return cluster, mult
    else:
        return None


def check_fa_rcm_consistency(fa_filename, rcm_filename):
    cluster_mult_rcm = defaultdict(int)
    num_rcm_reads = 0
    with open(rcm_filename) as rcm:
        for line in rcm:
            cluster = int(line.split("\t")[1])
            cluster_mult_rcm[cluster] += 1
            num_rcm_reads += 1

    num_fa_reads = 0
    is_ok = True
    with smart_open(fa_filename) as fa:
        for record in SeqIO.parse(fa, "fasta"):
            id = str(record.description)
            cluster, mult = parse_cluster_mult(id)
            if not cluster_mult_rcm[cluster] == mult:
                print id, cluster_mult_rcm[cluster]
                is_ok = False
            num_fa_reads += mult

    if not num_rcm_reads == num_fa_reads:
        is_ok = False
        print "Sizes are inconsistent: %d input reads and total multiplicity is %d" % (num_rcm_reads, num_fa_reads)

    return is_ok


def parse_cmd_line():
    parser = ArgumentParser(description="Compress equal clusters")
    parser.add_argument("input",
                        type=str,
                        help="input FASTA file")
    parser.add_argument("output",
                        type=str,
                        help="output FASTA file")
    parser.add_argument("--tmp-fa-file", "-T",
                        type=str,
                        default="",
                        help="temporary file for ig_trie_compressor output (default: <empty>)")
    parser.add_argument("--tmp-rcm-file", "-m",
                        type=str,
                        default="",
                        help="map file (default: <empty>)")
    parser.add_argument("--rcm", "-r",
                        type=str,
                        default="",
                        help="rcm file mapping reads out of scope to clusters from input file, leave empty to skip fixup stage (default: <empty>)")
    parser.add_argument("--output-rcm", "-R",
                        type=str,
                        default="",
                        help="output rcm file to combine --rcm file with --tmp-rcm-file, leave empty to rewrite existing --rcm file (default: <empty>)")

    args = parser.parse_args()

    if (not args.tmp_fa_file):
        args.tmp_fa_file = tempfile.mkstemp(suffix=".fa", prefix="igrec_")[1]
    if (not args.tmp_rcm_file):
        args.tmp_rcm_file = tempfile.mkstemp(suffix=".rcm", prefix="igrec_")[1]

    return args


if __name__ == "__main__":
    args = parse_cmd_line()

    print "Command line: %s" % " ".join(sys.argv)

    home_directory = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    run_trie_compressor = os.path.join(home_directory, 'build/release/bin/ig_trie_compressor')
    cmd_line = "%s -i %s -o %s -m %s -s true" % (run_trie_compressor, args.input, args.tmp_fa_file, args.tmp_rcm_file)
    os.system(cmd_line)

    input_read_id2mult = {}
    input_read_ids = []
    with open(args.input) as fin:
        for record in SeqIO.parse(fin, "fasta"):
            input_read_ids.append(record.description)
            cluster, mult = parse_cluster_mult(record.description)
            input_read_id2mult[record.description] = mult

    input_read_id2compressed_cluster_idx = rcm2rcmint(read_rcm_map(args.tmp_rcm_file))

    compressed_cluster_idx2mult = defaultdict(int)
    for input_read_id, compressed_cluster_idx in input_read_id2compressed_cluster_idx.iteritems():
        compressed_cluster_idx2mult[compressed_cluster_idx] += input_read_id2mult[input_read_id]

    # Fix IDs
    with smart_open(args.tmp_fa_file, "r") as fin, smart_open(args.output, "w") as fout:
        for i, record in enumerate(SeqIO.parse(fin, "fasta")):
            record.id = record.description = "cluster___%d___size___%d" % (i, compressed_cluster_idx2mult[i])
            SeqIO.write(record, fout, "fasta")

    # Fix RCM if provided
    if args.rcm:
        print "Fixing RCM file..."
        print "Checking input consistency..."
        if not check_fa_rcm_consistency(args.input, args.rcm):
            exit(1)

        if not args.output_rcm:
            args.output_rcm = args.rcm

        rcm = rcm2rcmint(read_rcm_list(args.rcm))

        with open(args.output_rcm, "w") as fout:
            for (outer_read_id, final_cluster_idx) in combine_rcms(rcm, input_read_ids, input_read_id2compressed_cluster_idx):
                fout.writelines("%s\t%s\n" % (outer_read_id, final_cluster_idx))

        print "Checking output consistency..."
        if not check_fa_rcm_consistency(args.output, args.output_rcm):
            exit(1)

    print "Removing temporary files: %s %s" % (args.tmp_fa_file, args.tmp_rcm_file)
    os.remove(args.tmp_fa_file)
    os.remove(args.tmp_rcm_file)
