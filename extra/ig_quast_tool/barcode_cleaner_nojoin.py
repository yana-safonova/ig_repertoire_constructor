#!/usr/bin/env python2

from argparse import ArgumentParser
from ig_basic import hamming_graph
from ig_basic import extract_barcode
from Bio import SeqIO
from collections import defaultdict


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
    g = hamming_graph(barcodes, tau=tau, method="knuth")
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





if __name__ == "__main__":
    parser = ArgumentParser(description="Fix barcode errors in dataset without joining; fixed protocol after discussion with Chudakov's team")
    parser.add_argument("input",
                        type=str,
                        help="input FASTQ file")
    parser.add_argument("output",
                        type=str,
                        help="output FASTQ file")
    parser.add_argument("--rcm", "-r",
                        type=str,
                        help="RCM file for final read-barcode map")
    parser.add_argument("--tau",
                        type=int,
                        default=3,
                        help="Minimum allowed distance between two barcodes [default %(default)s]")
    parser.add_argument("--no-fix-spaces", "-S",
                        action='store_true',
                        help="Do not replace spaces by underlines in read IDs")
    parser.add_argument("--min-size", "-m",
                        type=int,
                        default=0,
                        help="minimal barcode size [default %(default)d]")
    parser.add_argument("--bad-output", "-b",
                        type=str,
                        help="FASTQ file for bad-barcoded reads")
    # parser.add_argument("--subs-map", "-M",
    #                     type=str,
    #                     help="file for subs table")

    args = parser.parse_args()

    barcodes_count = defaultdict(int)

    print("Reading FASTQ...")
    with smart_open(args.input, "r") as fh:
        data = list(SeqIO.parse(fh, "fastq"))

    if not args.no_fix_spaces:
        for record in data:
            record.description = str(record.description).replace(" ", "_")
            record.id = record.name = record.description

    for record in data:
        barcode = extract_barcode(record.id)
        barcodes_count[barcode] += 1

    print("%d barcodes readed, %d unique barcodes" % (sum(barcodes_count.itervalues()),
                                                      len(barcodes_count)))

    barcodes = list(barcodes_count.keys())

    discarded_barcodes = set()

    print "Discarding rare barcodes for minimal size %d..." % args.min_size
    barcodes = discard_rare_barcodes(barcodes, args.min_size, barcodes_count, discarded_barcodes)
    print "Barcodes discarded: %d" % len(discarded_barcodes)

    for tau in range(1, args.tau):
        print "Discarding close barcodes for distance %d" % tau
        barcodes = discard_close_barcodes(barcodes, tau, barcodes_count, discarded_barcodes)
        print "Barcodes discarded: %d" % len(discarded_barcodes)


    print "Collecting output..."
    good_out = []
    bad_out = []
    for record in data:
        barcode = extract_barcode(record.id)
        if barcode not in discarded_barcodes:
            good_out.append(record)
        else:
            bad_out.append(record)

    print "Good reads: %d; bad reads: %d" % (len(good_out), len(bad_out))

    print "Writing output..."
    with smart_open(args.output, "w") as fh:
        SeqIO.write(good_out, fh, "fastq")

    if args.bad_output:
        with smart_open(args.bad_output, "w") as fh:
            SeqIO.write(bad_out, fh, "fastq")

    if args.rcm:
        with smart_open(args.rcm, "w") as fh:
            for record in good_out:
                barcode = extract_barcode(record.id)
                fh.write("%s\t%s\n" % (record.id, barcode))

    print "Done!"
