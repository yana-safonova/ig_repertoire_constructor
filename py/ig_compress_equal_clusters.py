#!/usr/bin/env python2

from argparse import ArgumentParser
from Bio import SeqIO
import sys
import os
import os.path
import tempfile
from collections import defaultdict

from ig_report_supernodes import smart_open


def parse_cluster_mult(id):
    import re
    id = str(id)
    m = re.match(r"^(intermediate_)?cluster___([0-9a-zA-Z]+)(___size___\d+)?___size___(\d+)$", id)
    if m:
        g = m.groups()
        cluster = g[1].strip()
        mult = int(g[3])
        return cluster, mult
    else:
        return None


assert parse_cluster_mult("cluster___10___size___20") == ('10', 20)
assert parse_cluster_mult("intermediate_cluster___10___size___20") == ('10', 20)
assert parse_cluster_mult("intermediate_cluster___10___size___5___size___20") == ('10', 20)
assert parse_cluster_mult("intermediatecluster___10___size___20") == None


def check_fa_rcm_consistency(fa_filename, rcm_filename):
    cluster_mult_rcm = defaultdict(int)
    num_rcm_reads = 0
    with smart_open(rcm_filename) as rcm:
        for line in rcm:
            cluster = line.split("\t")[1].strip()
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


if __name__ == "__main__":
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
    parser.add_argument("--tmp-map-file", "-m",
                        type=str,
                        default="",
                        help="map file (default: <empty>)")
    parser.add_argument("--rcm", "-r",
                        type=str,
                        default="",
                        help="rcm file for fixup, empty for skip this stage (default: <empty>)")
    parser.add_argument("--output-rcm", "-R",
                        type=str,
                        default="",
                        help="output rcm file for fixup, empty for rewriting existance file (default: <empty>)")
    parser.add_argument("--barcode-coverage-stats", "-S",
                        type=str,
                        default="",
                        help="file for output barcode coverage statistics; empty for non-producing (default: <empty>)")
    parser.add_argument("--barcode-threshold", "-b",
                        type=int,
                        default=1,
                        help="minimal number of different supporting barcodes (default: %(default)d)")
    parser.add_argument("--barcode", "-B",
                        action="store_true",
                        dest="count_barcodes",
                        help="add the number of different barcodes forming a cluster into cluster ID")
    parser.add_argument("--no-barcode",
                        action="store_false",
                        dest="count_barcodes",
                        help="do not add the number of different barcodes forming a cluster into cluster ID (default)")
    parser.set_defaults(count_barcodes=False)

    args = parser.parse_args()
    # Fix RCM if provided
    if args.rcm:
        print "Check input consistency..."
        if not check_fa_rcm_consistency(args.input, args.rcm):
            exit(1)

    if (not args.tmp_fa_file):
        args.tmp_fa_file = tempfile.mkstemp(suffix=".fa", prefix="igrec_")[1]
    if (not args.tmp_map_file):
        args.tmp_map_file = tempfile.mkstemp(suffix=".map", prefix="igrec_")[1]

    print "Command line: %s" % " ".join(sys.argv)

    home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/../'
    run_trie_compressor = os.path.join(home_directory, 'build/release/bin/./ig_trie_compressor')
    cmd_line = "%s -i %s -o %s -m %s" % (run_trie_compressor, args.input, args.tmp_fa_file, args.tmp_map_file)
    os.system(cmd_line)

    with smart_open(args.tmp_map_file) as fin:
        input_read_num2compressed_cluster = [line.strip() for line in fin]

    input_read_num2mult = []
    cluster2input_read_num = {}
    with smart_open(args.input) as fin:
        for i, record in enumerate(SeqIO.parse(fin, "fasta")):
            cluster, mult = parse_cluster_mult(record.description)
            cluster2input_read_num[cluster] = i
            input_read_num2mult.append(mult)

    compressed_cluster2mult = defaultdict(int)
    compressed_cluster2barcodes_number = defaultdict(int)
    for compressed_cluster, mult in zip(input_read_num2compressed_cluster, input_read_num2mult):
        compressed_cluster2mult[compressed_cluster] += mult
        compressed_cluster2barcodes_number[compressed_cluster] += 1

    if args.barcode_coverage_stats:
        with smart_open(args.tmp_fa_file, "r") as fin, smart_open(args.barcode_coverage_stats, "w") as fout:
            fout.write("%s\t%s\t%s\n" % ("cluster", "barcodes", "size"))
            for i, record in enumerate(SeqIO.parse(fin, "fasta")):
                n_barcodes = compressed_cluster2barcodes_number[str(i)]
                mult = compressed_cluster2mult[str(i)]
                fout.write("%d\t%d\t%d\n" % (i, n_barcodes, mult))

    # Fix IDs
    with smart_open(args.tmp_fa_file, "r") as fin, smart_open(args.output, "w") as fout:
        for i, record in enumerate(SeqIO.parse(fin, "fasta")):
            n_barcodes = compressed_cluster2barcodes_number[str(i)]
            if args.barcode_threshold > n_barcodes:
                continue
            if args.count_barcodes:
                record.id = record.description = "cluster___%dB%d___size___%d" % (i, n_barcodes, compressed_cluster2mult[str(i)])
            else:
                record.id = record.description = "cluster___%d___size___%d" % (i, compressed_cluster2mult[str(i)])
            SeqIO.write(record, fout, "fasta")

    # Fix RCM if provided
    if args.rcm:
        print "Fixing RCM file..."

        if not args.output_rcm:
            args.output_rcm = args.rcm

        with smart_open(args.rcm, "r") as fin:
            rcm = [line.strip().split("\t") for line in fin]
            rcm = [(id, cluster) for id, cluster in rcm]

        with smart_open(args.output_rcm, "w") as fout:
            for id, cluster in rcm:
                compressed_cluster = input_read_num2compressed_cluster[cluster2input_read_num[cluster]]
                n_barcodes = compressed_cluster2barcodes_number[compressed_cluster]
                if args.barcode_threshold > n_barcodes:
                    continue
                if args.count_barcodes:
                    compressed_cluster += "B" + str(n_barcodes)
                fout.writelines("%s\t%s\n" % (id, compressed_cluster))

        print "Check output consistency..."
        if not check_fa_rcm_consistency(args.output, args.output_rcm):
            exit(1)

    print "Remove temporary files: %s %s" % (args.tmp_fa_file, args.tmp_map_file)
    os.remove(args.tmp_fa_file)
    os.remove(args.tmp_map_file)
