#!/usr/bin/env python2

from argparse import ArgumentParser
from Bio import SeqIO
import sys


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


def parse_size(s):
    import re

    m = re.match(r".*___size___(\d+)", s)

    if m:
        g = m.groups()
        return int(g[0])
    else:
        return 1

# assert(parse_abundance("dsdsfsd_abundance:10") == 10)
assert(parse_size("dsdsfsd___size___10") == 10)

def parse_multiplicity(s):
    import re

    m = re.match(r".*_multiplicity_(\d+)_", s)

    if m:
        g = m.groups()
        return int(g[0])
    else:
        return 1

assert parse_multiplicity("antibody_1_multiplicity_999_copy_1") == 999



if __name__ == "__main__":
    parser = ArgumentParser(description="Report supernodes")
    parser.add_argument("input",
                        type=str,
                        help="input FASTA file with abundances in ids")
    parser.add_argument("output",
                        type=str,
                        help="output FASTA file")
    parser.add_argument("--limit", "-l",
                        type=int,
                        default=5,
                        help="abundance limit (default: %(default)s)")

    args = parser.parse_args()

    print "Supernode reporter started..."
    print "Command line: %s" % " ".join(sys.argv)

    result = []
    with smart_open(args.input, "r") as fin:
        for record in SeqIO.parse(fin, "fasta"):
            abundance = parse_size(str(record.id))
            if abundance >= args.limit:
                old_id = str(record.id)
                if "__size__" in old_id:
                    new_id = old_id
                else:
                    new_id = "%s___size___%d" % (old_id, abundance)
                record.id = record.name = record.description = new_id
                result.append(record)

    with smart_open(args.output, "w") as fout:
        SeqIO.write(result, fout, "fasta")

    print "Supernode reporter done"
