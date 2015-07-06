#!/usr/bin/env python2

from argparse import ArgumentParser
from ig_basic import hamming_graph
# import glob
# import os
from Bio import SeqIO
from read_barcode_splitter import extract_barcode
from collections import defaultdict


# import resource, sys
# resource.setrlimit(resource.RLIMIT_STACK, (2**29,-1))
# sys.setrecursionlimit(10**8)


import contextlib
@contextlib.contextmanager
def smart_open(filename, mode="r"):
    """
    From http://stackoverflow.com/questions/17602878/how-to-handle-both-with-open-and-sys-stdout-nicely
    """
    import gzip
    import re
    from sys import stdout, stdin

    if "w" in mode:
        MODE = "w"
    elif "a" in mode:
        MODE= "a"
    else:
        MODE = "r"

    if filename != '-':
        if re.match(r"^.*\.gz$", filename):
            assert(MODE != "a")
            fh = gzip.open(filename, mode=MODE)
        else:
            fh = open(filename, mode=mode)
    else:
        assert(MODE != "a")
        fh = stdout if MODE == "w" else stdin
    try:
        yield fh
    finally:
        if fh is not stdout and fh is not stdin:
            fh.close()


if __name__ == "__main__":
    parser = ArgumentParser(description="Fix barcode errors in dataset")
    parser.add_argument("input",
                        type=str,
                        help="input FASTQ file")
    parser.add_argument("output",
                        type=str,
                        help="output FASTQ file")
    barcode_list_group = parser.add_mutually_exclusive_group(required=True)
    barcode_list_group.add_argument("--minimal-size", "-m",
                                    type=int,
                                    help="minimal barcode size")
    barcode_list_group.add_argument("--barcode-list", "-b",
                                    type=str,
                                    help="original barcode list")
    parser.add_argument("--tau",
                        type=int,
                        default=3,
                        help="Maximum corrected error [default %(default)s]")
    parser.add_argument("--bad-output", "-B",
                        type=str,
                        help="FASTQ file for bad-barcoded reads")
    parser.add_argument("--supernodes", "-S",
                        type=str,
                        help="file for estimated original barcode list")
    parser.add_argument("--hgraph-original",
                        type=str,
                        help="file for original barcodes H-graph (in .dot format)")
    parser.add_argument("--hgraph-data",
                        type=str,
                        help="file for data barcodes H-graph (in .dot format)")
    parser.add_argument("--subs-map", "-M",
                        type=str,
                        help="file for subs table")

    args = parser.parse_args()

    barcodes_count = defaultdict(int)

    print("Reading FASTQ...")
    with smart_open(args.input, "r") as fh:
        for record in SeqIO.parse(fh, "fastq"):
            barcode = extract_barcode(record.id)
            barcodes_count[barcode] += 1

    print("%d barcodes readed, %d unique barcodes" % (sum(barcodes_count.itervalues()),
                                                      len(barcodes_count)))

    if args.minimal_size is not None:
        original_barcodes = [barcode for barcode, n in barcodes_count.iteritems() if n >= args.minimal_size]
        if args.supernodes is not None:
            with smart_open(args.supernodes, "w") as fh:
                fh.writelines([barcode + "\n" for barcode in original_barcodes])
    elif args.barcode_list is not None:
        with smart_open(args.barcode_list, "r") as fh:
            original_barcodes = [barcode.strip() for barcode in fh]

    original_barcodes = set(original_barcodes)
    data_barcodes = set(barcodes_count.iterkeys())

    print "Original barcodes %d" % len(original_barcodes)
    print "Data barcodes %d" % len(data_barcodes)

    from knuth_index import KnuthIndex
    from Levenshtein import hamming

    original_barcode_mult = [barcodes_count[barcode] if barcode in barcodes_count else 0 for barcode in original_barcodes]

    print "Knuth k-mer index construction..."
    tree = KnuthIndex(list(original_barcodes),
                      tau=args.tau,
                      priority=original_barcode_mult)
    print "Index constructed"

    bad_barcodes = []
    barcode_barcode = {}
    dists = []

    for barcode in data_barcodes:
        if barcode in original_barcodes:
            barcode_barcode[barcode] = barcode
            dists.append(0)
        else:
            neib = tree.get_nn(barcode, nonn="None")
            if neib is not None:
                barcode_barcode[barcode] = neib
                dists.append(hamming(neib, barcode))
            else:
                bad_barcodes.append(barcode)
                dists.append(-1)

    print "Bad-coded reads %d unique barcodes %d" % (sum(barcodes_count[barcode] for barcode in bad_barcodes),
                                                     len(bad_barcodes))

    print "Well-coded reads %d unique barcodes %d" % (sum(barcodes_count[barcode] for barcode in barcode_barcode.iterkeys()),
                                                      len(barcode_barcode))


    dist_hist = defaultdict(int)
    dist_barcodes = defaultdict(list)
    for dist, barcode in zip(dists, data_barcodes):
        dist_hist[dist] += 1
        dist_barcodes[dist].append(barcode)

    for tau in range(-1, args.tau + 1):
        print "Barcodes with dist %d: %d, reads %d" % (tau,
                                                       dist_hist[tau],
                                                       sum(barcodes_count[barcode] for barcode in dist_barcodes[tau]))

    print("Reading and writing FASTQ...")
    output = []
    bad_output = []

    # TODO Do it more accurately
    def change_barcode(s, old_barcode, new_barcode):
        return s.replace("UMI:" + old_barcode + ":",
                         "UMI:" + new_barcode + ":",
                         1)

    with smart_open(args.input, "r") as fh:
        for record in SeqIO.parse(fh, "fastq"):
            barcode = extract_barcode(record.id)
            if barcode in barcode_barcode:
                record.id = change_barcode(record.id, barcode, barcode_barcode[barcode])
                output.append(record)
            else:
                bad_output.append(record)

    with smart_open(args.output, "w") as fh:
        SeqIO.write(output, fh, "fastq")

    if args.bad_output is not None:
        with smart_open(args.bad_output, "w") as fh:
            SeqIO.write(bad_output, fh, "fastq")

    if args.hgraph_data is not None:
        print "Construction H-graph for data barcodes..."
        g_data = hamming_graph(list(data_barcodes), tau=args.tau)
        g_data.write_dot(args.hgraph_data)
        print "Constructed"

    if args.hgraph_original is not None:
        print "Construction H-graph for original barcodes..."
        g_original = hamming_graph(list(original_barcodes), tau=args.tau)
        from scipy.stats import itemfreq
        if len(g_original.es) > 0:
            print itemfreq(g_original.es["weight"])
        g_original.write_dot(args.hgraph_original)
        print "Constructed"

    if args.subs_map is not None:
        with smart_open(args.subs_map, "w") as fh:
            for _f, _t in barcode_barcode.iteritems():
                fh.write("%s %s\n" % (_f, _t))
