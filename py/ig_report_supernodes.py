#!/usr/bin/env python2

from argparse import ArgumentParser
from Bio import SeqIO
import sys


import os
current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = current_dir + "/../"
sys.path.append(igrec_dir + "/py/utils")
sys.path.append(igrec_dir + "/py/pipeline/")
import support
sys.path.append(igrec_dir + "/py/")
from ash_python_utils import idFormatByFileName, smart_open

def parse_size(s):
    import re

    m = re.match(r"^.*___size___(\d+)$", s)

    if m:
        g = m.groups()
        return int(g[0])
    else:
        return None

assert(parse_size("dsdsfsd___size___10") == 10)


if __name__ == "__main__":
    parser = ArgumentParser(description="Report supernodes")
    parser.add_argument("input",
                        type=str,
                        help="input FASTA/FASTQ file with abundances in ids")
    parser.add_argument("output",
                        type=str,
                        help="output FASTA/FASTQ file")
    parser.add_argument("--limit", "-l",
                        type=int,
                        default=5,
                        help="size limit (default: %(default)s)")

    args = parser.parse_args()

    print "Supernode reporter started..."
    print "Command line: %s" % " ".join(sys.argv)

    input_size = output_size = 0
    with smart_open(args.input, "r") as fin, smart_open(args.output, "w") as fout:
        for record in SeqIO.parse(fin, idFormatByFileName(args.input)):
            input_size += 1
            id = str(record.description)
            size = parse_size(id)
            assert id is not None
            if size >= args.limit:
                SeqIO.write(record, fout, idFormatByFileName(args.output))
                output_size += 1


    print "%d antibody clusters have abundance >= %d" % (output_size, args.limit)
    print "%d lowly abundant antibody clusters will be discarded" % (input_size - output_size, )
    print "Highly abundant clusters were written to " + args.output
    print "Supernode reporter done"
