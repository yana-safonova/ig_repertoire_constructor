#!/usr/bin/env python2

from simulate import run_mixcr
from argparse import ArgumentParser

if __name__ == "__main__":
    parser = ArgumentParser(description="MiXCR wrapper")
    parser.add_argument("input",
                        type=str,
                        help="input FASTA/FASTQ file with MiSeq Rep-Seq single reads")
    parser.add_argument("output",
                        type=str,
                        help="output dir")
    parser.add_argument("--loci", "-l",
                        type=str,
                        default="IGH",
                        help="loci: %(default)s)")

    args = parser.parse_args()

    run_mixcr(args.input, args.output, loci=args.loci)
