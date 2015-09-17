#!/usr/bin/env python2

from argparse import ArgumentParser
from barcode_cleaner import smart_open

from ig_basic import extract_barcode, most_popular_element
from collections import defaultdict
import numpy as np


def consensus(reads):
    import numpy as np

    n = min(map(len, reads))

    reads = reads[:]
    reads = [read for read in reads if len(read) > 0]
    if len(reads) == 0 or n == 0:
        return None

    for i in xrange(len(reads)):
        reads[i] = reads[i][:n]

    nucleotides = np.array(['A', 'C', 'G', 'T'])
    nuc_ind = {nucleotides[i]: i for i in range(len(nucleotides))}

    consm = np.zeros((4, n), dtype=int)
    for j in xrange(len(reads)):
        s = reads[j]
        for i in xrange(n):
            consm[nuc_ind[s[i]], i] += 1

    consensus = "".join(nucleotides[np.argmax(consm, axis=0)])

    return consensus


def print_barcode_abundance(barcode2abundance, barcodes):
    barcode_abundances = [barcode2abundance[_] for _ in barcodes]

    print "min = %d\t\tmax = %d\tmean = %0.2f" % (min(barcode_abundances),
                                           max(barcode_abundances),
                                           np.mean(barcode_abundances))

if __name__ == "__main__":
    parser = ArgumentParser(description="Compare rcm-file")
    parser.add_argument("--input", "-i",
                        type=str,
                        required=True,
                        help="Input RCM file")
    parser.add_argument("--barcode_rcm", "-b",
                        type=str,
                        help="Barcode RCM file")
    # parser.add_argument("-reads", "-s",
    #                     type=str,
    #                     required=False,
    #                     help="Input (cropped) reads file (FASTA)")
    parser.add_argument("--rate", "-r",
                        type=float,
                        default=0.5,
                        help="<<good>> barcode coverage rate threshold (default = %(default)0.2f})")

    args = parser.parse_args()
    ids = []
    cliques = []
    barcodes = []

    with smart_open(args.input) as fh:
        for l in fh:
            id, clique = l.split("\t")
            clique = int(clique)
            ids.append(id)
            cliques.append(clique)

    if args.barcode_rcm is None:
        for id in ids:
            barcodes.append(extract_barcode(id))
    else:
        id2barcode = {}
        with smart_open(args.barcode_rcm) as barcode_rcm:
            for l in barcode_rcm:
                id, barcode = l.split("\t")
                id = id.strip()
                barcode = barcode.strip()
                id2barcode[id] = barcode
        barcodes = [id2barcode[_id] for _id in ids]

    # reads = []
    # read_ids = []
    # id2i = {}
    #
    # print "Reading input reads..."
    # with smart_open(args.reads, "r") as fh:
    #     for record in SeqIO.parse(fh, "fasta"):
    #         id = str(record.id)
    #         read = str(record.seq)
    #         reads.append(read)
    #         read_ids.append(id)
    #
    # for id, id_read in zip(ids, read_ids):
    #     assert(id == id_read)
    #
    #
    # for i in xrange(len(ids)):
    #     id2i[ids[i]] = i

    id2clique = {id: clique for id, clique in zip(ids, cliques)}

    barcode2cliques = defaultdict(list)
    for barcode, clique in zip(barcodes, cliques):
        barcode2cliques[barcode].append(clique)

    barcode2ids = defaultdict(list)
    for barcode, id in zip(barcodes, ids):
        barcode2ids[barcode].append(id)

    # barcode2consensus = {}
    # barcode2dists = {}
    # for barcode, ii in barcode2ids.iteritems():
    #     local_reads = [reads[id2i[i]] for i in ii]
    #     l = min(map(len, local_reads))
    #     trimmed_reads = [read[:l] for read in local_reads]
    #     barcode2consensus[barcode] = cons = consensus(trimmed_reads)
    #     barcode2dists[barcode] = [levenshtein(read, cons) for read in trimmed_reads]
    #
    # for barcode, dists in barcode2dists.iteritems():
    #     print barcode, dists

    barcode2mp_clique = {id: most_popular_element(cliques) for id, cliques in barcode2cliques.iteritems()}
    barcode2abundance = {id: len(cliques) for id, cliques in barcode2cliques.iteritems()}

    good_barcodes = []
    bad_barcodes = []
    for barcode, cliques in barcode2cliques.iteritems():
        mp_clique = most_popular_element(cliques)
        nmp = sum([clique == mp_clique for clique in cliques])

        rate = float(nmp) / len(cliques)
        if rate >= args.rate:
            good_barcodes.append(barcode)
        else:
            bad_barcodes.append(barcode)

    print "Good barcodes: %d" % len(good_barcodes)
    print_barcode_abundance(barcode2abundance, good_barcodes)
    print "Bad barcodes: %d" % len(bad_barcodes)
    print_barcode_abundance(barcode2abundance, bad_barcodes)

    # sys.exit()
    #
    #
    # print "Reading input reads..."
    # with smart_open(args.reads, "r") as fh:
    #     id2read = {str(record.id): str(record.seq)
    #                for record in SeqIO.parse(fh, "fasta")}
    #
    # barcode2ids = defaultdict(list)
    # for id in id2clique.keys():
    #     barcode = extract_barcode(id)
    #     barcode2ids[barcode].append(id)
    #
    # barcode2consensus = {}
    # barcode2dists = {}
    # for barcode, ids in barcode2ids.iteritems():
    #     cons = barcode2consensus[barcode] = consensus([id2read[id] for id in ids])
    #     dists = [levenshtein(cons, id2read[id]) for id in ids]
    #     barcode2dists[barcode] = dists
    #
    # for barcode in good_barcodes:
    #     print barcode2dists[barcode]
    # for barcode in bad_barcodes:
    #     print barcode2dists[barcode]
