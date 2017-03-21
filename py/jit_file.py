#!/usr/bin/env python2

from argparse import ArgumentParser
from simulate import jit_fx_file
import sys

if __name__ == "__main__":
    parser = ArgumentParser(description="Add random sequencing errors to FASTA/FASTQ file")
    parser.add_argument("input",
                        type=str,
                        help="input FASTA/FASTQ file")
    parser.add_argument("output",
                        type=str,
                        help="output FASTA/FASTQ file")
    parser.add_argument("--error-rate",
                        type=float,
                        default=0,
                        help="error rate (default: %(default)0.2f)")
    parser.add_argument("--seed", "-S",
                        type=int,
                        default=0,
                        help="random seed (default: %(default)d)")

    args = parser.parse_args()
    print "FASTA/FASTQ file jittering"
    print "Command line: %s" % " ".join(sys.argv)

    jit_fx_file(args.input, args.output, error_rate=args.error_rate,
                random_errors=True,
                min_error=0,
                erroneous_site_len=10005000, seed=args.seed)
    print "Jittering done!"
