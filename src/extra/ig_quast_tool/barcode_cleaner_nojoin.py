#!/usr/bin/env python2

from argparse import ArgumentParser
from ig_basic import hamming_graph
from ig_basic import extract_barcode
from Bio import SeqIO
from collections import defaultdict

import os
import os.path
import sys
import matplotlib
matplotlib.use('Agg')
import seaborn as sns

current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = current_dir + "/../../../"
sys.path.append(igrec_dir + "/py/utils")
sys.path.append(igrec_dir + "/py/pipeline/")
sys.path.append(igrec_dir + "/py/")
from ash_python_utils import idFormatByFileName, smart_open

sys.path.append(igrec_dir + "/py")
from igquast_impl import consensus
from igquast_impl import initialize_plot, save_plot


def discard_rare_barcodes(barcodes, min_size, barcode2size, discarded):
    result = []
    for barcode in barcodes:
        if barcode2size[barcode] < min_size:
            discarded.add(barcode)
        else:
            result.append(barcode)

    return result


def discard_close_barcodes(barcodes, tau, barcode2size, discarded):
    result = []
    g = hamming_graph(barcodes, tau=tau, method="bk")
    vs = list(g.vs())
    vs.sort(key=lambda v: barcode2size[v["read"]])
    for v in vs:
        if v["read"] in discarded:
            next
        discard = False
        for vv in v.neighbors():
            if vv["read"] not in discarded:
                discard = True
                discarded.add(v["read"])
                break
        if not discard:
            result.append(v["read"])

    return result


def hamming(X, Y):
    res = 0
    for x, y in zip(X, Y):
        if x != y:
            res += 1
    return res


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Fix barcode errors in dataset without joining; fixed protocol after discussion with Chudakov's team")
    parser.add_argument("input", type=str, help="input FASTA/FASTQ file")
    parser.add_argument("output", type=str, help="output FASTA/FASTQ file")
    parser.add_argument("--rcm",
                        "-r",
                        type=str,
                        help="RCM file for final read-barcode map")
    parser.add_argument("--tau",
                        type=int,
                        default=3,
                        help="Minimum allowed distance between two barcodes [default %(default)s]")
    parser.add_argument("--no-fix-spaces",
                        "-S",
                        action='store_true',
                        help="Do not replace spaces by underlines in read IDs")
    parser.add_argument("--min-size",
                        "-m",
                        type=int,
                        default=0,
                        help="minimal barcode size [default %(default)d]")
    parser.add_argument("--bad-output",
                        "-b",
                        type=str,
                        help="FASTA/FASTQ file for bad-barcoded reads")
    parser.add_argument("--distance-plot",
                        type=str,
                        default="",
                        help="Figure file for distance distribution plot")
    parser.add_argument("--distance-threshold",
                        "-d",
                        type=int,
                        default=10,
                        help="distance threshold [default %(default)d]")
    parser.add_argument("--lengths",
                        type=str,
                        help="file for read length stats")
    # parser.add_argument("--subs-map", "-M",
    #                     type=str,
    #                     help="file for subs table")

    args = parser.parse_args()

    barcodes_count = defaultdict(int)

    print("Reading library...")
    with smart_open(args.input, "r") as fh:
        data = list(SeqIO.parse(fh, idFormatByFileName(args.input)))

    if not args.no_fix_spaces:
        for record in data:
            record.description = str(record.description).replace(" ", "_")
            record.id = record.name = record.description

    # Omit reads with Ns
    data = [record for record in data if record.seq.count("N") == 0]

    data = [record for record in data if extract_barcode(record.id) is not None]

    clusters = defaultdict(list)

    for record in data:
        barcode = extract_barcode(record.id)
        barcodes_count[barcode] += 1
        clusters[barcode].append(record)

    cluster_consensus = {}
    distances = []
    print "Compute consensus by barcodes"
    for barcode, reads in clusters.iteritems():
        reads = [str(read.seq) for read in reads]
        cons = cluster_consensus[barcode] = consensus(reads)
        distances += [hamming(cons, read) for read in reads]

    if args.distance_plot:
        initialize_plot()
        sns.distplot(distances)
        save_plot(args.distance_plot, ("png", "pdf"))
        print "distance plot saved to %s.{png,pdf}" % args.distance_plot

    print("%d barcodes readed, %d unique barcodes" %
          (sum(barcodes_count.itervalues()), len(barcodes_count)))

    barcodes = list(barcodes_count.keys())

    discarded_barcodes = set()

    print "Discarding rare barcodes for minimal size %d..." % args.min_size
    barcodes = discard_rare_barcodes(barcodes, args.min_size, barcodes_count,
                                     discarded_barcodes)
    print "Barcodes discarded: %d" % len(discarded_barcodes)

    if args.tau:
        tau = args.tau - 1
        print "Discarding close barcodes for distance %d" % tau
        barcodes = discard_close_barcodes(barcodes, tau, barcodes_count,
                                          discarded_barcodes)
        print "Barcodes discarded: %d" % len(discarded_barcodes)

    print "Collecting output..."
    good_out = []
    bad_out = []
    for record in data:
        barcode = extract_barcode(record.id)
        if barcode not in discarded_barcodes and hamming(
                record.seq,
                cluster_consensus[barcode]) <= args.distance_threshold:
            good_out.append(record)
        else:
            bad_out.append(record)

    print "Good reads: %d; bad reads: %d" % (len(good_out), len(bad_out))

    print "Writing output..."
    with smart_open(args.output, "w") as fh:
        SeqIO.write(good_out, fh, idFormatByFileName(args.output))

    if args.bad_output:
        with smart_open(args.bad_output, "w") as fh:
            SeqIO.write(bad_out, fh, idFormatByFileName(args.bad_output))

    if args.rcm:
        with smart_open(args.rcm, "w") as fh:
            for record in good_out:
                barcode = extract_barcode(record.id)
                fh.write("%s\t%s\n" % (record.id, barcode))

    barcode2lengths = defaultdict(list)
    if args.lengths:
        with smart_open(args.lengths, "w") as fh:
            for record in good_out:
                barcode = extract_barcode(record.id)
                l = len(record.seq)
                barcode2lengths[barcode].append(l)

            for barcode, lst in barcode2lengths.iteritems():
                count = defaultdict(int)
                for l in lst:
                    count[l] += 1
                if len(count) > 1:
                    items = count.items()
                    items.sort(key=lambda x: x[1], reverse=True)
                    assert items[0][1] >= items[1][1]
                    if items[1][1] >= 3:
                        out = ", ".join(["%d: %d" % (l, c) for l, c in items])
                        fh.write("%s: %s\n" % (barcode, out))
                        for read in clusters[barcode]:
                            SeqIO.write(read, fh, "fasta")
    print "Done with " + args.input
