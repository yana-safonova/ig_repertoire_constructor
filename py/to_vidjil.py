#!/usr/bin/env python2

import argparse
import os
import json
import csv
import sys
import pandas as pd
import numpy as np
import datetime
from Bio import SeqIO
from collections import defaultdict
import scipy
import scipy.stats

current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = current_dir + "/../"
sys.path.append(igrec_dir + "/py/")
from ash_python_utils import idFormatByFileName, smart_open
from ig_compress_equal_clusters import parse_cluster_mult
import build_info


def compute_average_read_length(initial_reads_file, rcm_file, fix_spaces=True, trim_trailing_underscores=True):
    # parse rcm and compute _average_read_length
    import re

    normalize_id = lambda s: s.replace(" ", "_") if fix_spaces else lambda s: s

    with smart_open(initial_reads_file) as fin:
        id2len = {normalize_id(str(record.description)): len(record.seq) for record in SeqIO.parse(fin, idFormatByFileName(args.initial_reads))}

    cluster2reads = defaultdict(list)
    with open(rcm_file) as rcm:
        for line in rcm:
            read_id, cluster_id = line.split()
            if trim_trailing_underscores:
                read_id = re.sub(r"_+$", "", read_id)
            cluster2reads[cluster_id].append(read_id)

    cluster2avlen = {cluster_id: np.mean([id2len[read_id] for read_id in read_ids]) for cluster_id, read_ids in cluster2reads.iteritems()}
    return cluster2avlen


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", help="input dir with IgReC results", required=True)
    parser.add_argument("-o", "--output", help="output .vidjil file", required=True)
    parser.add_argument("--initial-reads", help="initial sequencing reads file [default: <input>/merged_reads.fastq]", required=False)
    args = parser.parse_args()

    if args.initial_reads is None:
        args.initial_reads = args.input + "/merged_reads.fastq"

    cluster2avlen = compute_average_read_length(args.initial_reads, args.input + "/final_repertoire.rcm")

    loci = []
    with open(args.input + "/vj_finder/alignment_info.csv") as tsvin:
        for row in csv.reader(tsvin, delimiter='\t'):
            loci.append(row[1])

    segmented = len(loci)

    from collections import Counter
    germline = {k: [v] for k, v in Counter(loci).iteritems()}

    with open(args.input + "/vj_finder/filtering_info.csv") as tsvin:
        nonsegmented = sum(1 for _ in tsvin)

    total = segmented + nonsegmented

    cdr_details = pd.read_csv(args.input + "/divan/cdr_details.txt", sep='\t')
    cdr_details.set_index("Clone_name", inplace=True)

    clones = []
    with open(args.input + "/final_repertoire.fa") as fin:
        for record in SeqIO.parse(fin, "fasta"):
            id = str(record.description)
            try:
                row = cdr_details.loc[id]
            except:
                print "Cluster %s is not identified by DivAn" % id
                row = defaultdict(str)
            cluster_id, mult = parse_cluster_mult(record.description)
            name = row["AA_seq"]  # TODO USe CDR + V + J format like in MiXCR; it requires divan output table extension
            clone = {"id": id,
                     "name": name,
                     "sequence": str(record.seq),
                     "reads": [mult],
                     "_average_read_length": [cluster2avlen[cluster_id]],
                     "germline": row["Chain_type"],
                     "seg": {
                         "5": {"name": row["V_hit"]},
                         "3": {"name": row["J_hit"]},
                         "cdr3": {"start": row["CDR3_start"], "stop": row["CDR3_end"], "seq": row["CDR3_nucls"]}
                        }
                     }
            clones.append(clone)

    sizes = [clone["reads"][0] for clone in clones]
    ranks = scipy.stats.rankdata(-np.array(sizes), method="ordinal")
    for clone, rank in zip(clones, ranks):
        clone["top"] = int(rank)

    command_line = "N/A"
    with open(args.input + "/igrec.log") as fin:
        prefix = "Command line: "
        for line in fin:
            if line.startswith(prefix):
                command_line = line[len(prefix):].strip()
                break

    # TODO export args into json on igrec.py run

    output = {"producer": "IgReC %s, git hash %s" % (build_info.version, build_info.git_hash),
              "timestamp": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),  # TODO get timestamp from igrec.json
              "vidjil_json_version": "2016b",
              "commandline": [command_line],

              "samples": {
                  "number": 1,
                  "original_names": [args.initial_reads],
                },

              "reads": {
                  "total": [total],
                  "segmented": [segmented],
                  "germline": germline
                },

              "clones": clones
              }

    with open(args.output, "w") as fout:
        json.dump(output, fout, sort_keys=True, indent=4)
