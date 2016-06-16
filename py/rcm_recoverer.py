#!/usr/bin/env python2

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


def generate_rcm(compressed_rcm_file_name, clique_ids_file_name, out_file):
    with smart_open(clique_ids_file_name, "r") as fh:
        compressed_read2clique = [int(s) for s in fh]

    with smart_open(compressed_rcm_file_name, "r") as fh:
        rcm = [line.split("\t") for line in fh]
    ids, cliques = zip(*[(id, int(clique)) for id, clique in rcm])

    with smart_open(out_file, "w") as fh:
        for i in xrange(len(ids)):
            fh.write("%s\t%d\n" % (ids[i], compressed_read2clique[cliques[i]]))


if __name__ == "__main__":
    parser = ArgumentParser(description="Recover .RCM file from .map file")
    parser.add_argument("--compressed-map", "-c",
                        type=str,
                        required=True,
                        help="rcm output file of trie compressor")
    parser.add_argument("--cliques-map", "-q",
                        type=str,
                        required=True,
                        help="cliques map")
    parser.add_argument("--output", "-o",
                        type=str,
                        required=True,
                        help="output rcm file")

    args = parser.parse_args()

    generate_rcm(args.compressed_map, args.cliques_map, args.output)
