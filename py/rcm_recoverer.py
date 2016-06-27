#!/usr/bin/env python2

from Bio import SeqIO
from argparse import ArgumentParser


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


def generate_rcm(reads_file_name, compressed_file_name, cliques_ids_file_name, out_file):
    # Obtain read ids
    with smart_open(reads_file_name, "r") as fh:
        ids = [str(record.id) for record in SeqIO.parse(fh, "fasta")]

    # Obtain compread2clique
    with smart_open(cliques_ids_file_name, "r") as fh:
        compread2clique = [int(s) for s in fh]

    with smart_open(compressed_file_name, "r") as fh:
        idmap = [int(s) for s in fh]

    with smart_open(out_file, "w") as fh:
        for i in xrange(len(ids)):
            fh.write("%s\t%d\n" % (ids[i], compread2clique[idmap[i]]))


if __name__ == "__main__":
    parser = ArgumentParser(description="Recover .RCM file from .map file")
    parser.add_argument("--ids-file", "-i",
                        type=str,
                        required=True,
                        help="FASTA file with ids")
    parser.add_argument("--output", "-o",
                        type=str,
                        required=True,
                        help="output RCM file")
    parser.add_argument("--compressed-map", "-c",
                        type=str,
                        required=True,
                        help="compressed map")
    parser.add_argument("--cliques-map", "-q",
                        type=str,
                        required=True,
                        help="cliques map")

    args = parser.parse_args()

    generate_rcm(args.ids_file, args.compressed_map, args.cliques_map, args.output)
