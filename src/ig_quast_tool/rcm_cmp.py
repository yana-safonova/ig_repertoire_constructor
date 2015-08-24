#!/usr/bin/env python2

from Bio import SeqIO
from argparse import ArgumentParser
from barcode_cleaner import smart_open

from ig_basic import extract_barcode, most_popular_element
from collections import defaultdict
import numpy as np

if __name__ == "__main__":
    parser = ArgumentParser(description="Compare rcm-file")
    parser.add_argument("--input", "-i",
                        type=str,
                        required=True,
                        help="Input RCM file")
    parser.add_argument("-reads", "-s",
                        type=str,
                        required=True,
                        help="Input reads file (FASTA)")
    parser.add_argument("--rate", "-r",
                        type=float,
                        default=0.9,
                        help="<<good>> barcode coverage rate threshold (default = %{default})")

    args = parser.parse_args()
    with smart_open(args.input) as fh:
        id2clique = {id: int(clique) for id, clique in [l.split("\t") for l in fh]}

    barcode2cliques = defaultdict(list)
    for id, clique in id2clique.iteritems():
        barcode = extract_barcode(id)
        barcode2cliques[barcode].append(clique)


    for barcode, cliques in barcode2cliques.iteritems():
        print barcode, cliques

    barcode2mp_clique = {id: most_popular_element(cliques) for id, cliques in barcode2cliques.iteritems()}
    barcode2abundance = {id: len(cliques) for id, cliques in barcode2cliques.iteritems()}

    good_barcodes = []
    bad_barcodes = []
    for barcode, cliques in barcode2cliques.iteritems():
        mp_clique = most_popular_element(cliques)
        nmp = sum([clique == mp_clique for clique in cliques])

        rate = float(nmp) / len(clique)
        if rate >= args.rate:
            good_barcodes.append(barcode)
        else:
            bad_barcodes.append(barcode)

    print "Good barcodes: %d" % len(good_barcodes)
    print "Bad barcodes: %d" % len(bad_barcodes)

    bad_barcode_abundances = [barcode2abundance[_] for _ in bad_barcodes]
    print "Bad barcode abundance: min = %d, max = %d, mean=%f" % (min(bad_barcode_abundances),
                                                                  max(bad_barcode_abundances),
                                                                  np.mean(bad_barcode_abundances))

    print "Reading input reads..."
    with smart_open(args.reads, "r") as fh:
        id2read = {str(record.id): record for record in seqIO.parse(fh, "fasta")}








