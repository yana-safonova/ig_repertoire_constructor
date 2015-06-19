#!/usr/bin/env python2

# from Bio import SeqIO
# import argparse
import glob
import os
import re


if __name__ == "__main__":
    try:
        os.mkdir("s3_bcodes_align")
    except Exception:
        print("Directory has already existed")

    try:
        os.remove("consensus.fasta")
    except Exception:
        pass
    try:
        os.remove("stats.txt")
    except Exception:
        pass


    for fname in open("list.txt"):
        fname = fname.strip()
        fname_aligned = re.sub(r"s3_bcodes", "s3_bcodes_align", fname)
        os.system("ashtools/align.py -d 0 -v --stats stats.txt --consensus consensus.fasta %s %s" % (fname, fname_aligned))
