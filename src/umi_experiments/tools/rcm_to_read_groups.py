import os
import shutil
import sys
import strop
from collections import defaultdict
from Bio import SeqIO


def ParseCommandLineParams():
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Create separate fasta files with reads from each cluster.")
    parser.add_argument("-r", "--rcm", type=str, dest="rcm_path", help="Rcm file path")
    parser.add_argument("-i", "--input", type=str, dest="reads_path", help="Path to file with input reads")
    parser.add_argument("-o", "--output", type=str, dest="output_path", help="Output directory path")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    params = parser.parse_args()

    return params


def ReadRcm(rcm_path):
    with open(rcm_path, "r") as rcm_file:
        ids_and_clusters = [strop.strip(line).split('\t') for line in rcm_file]
        return zip(*ids_and_clusters)


def ReadFasta(reads_path):
    with open(reads_path, "r") as reads_file:
        return SeqIO.to_dict(SeqIO.parse(reads_file, "fasta"))


def main():
    params = ParseCommandLineParams()

    read_ids, cluster_idxs = ReadRcm(params.rcm_path)
    reads = ReadFasta(params.reads_path)
    print "Done reading"

    cluster_to_read_ids = defaultdict(list)
    for i in range(len(read_ids)):
        cluster_to_read_ids[cluster_idxs[i]].append(read_ids[i])
    print "Done restoring clusters"

    if os.path.exists(params.output_path):
        shutil.rmtree(params.output_path)
    os.makedirs(params.output_path)
    all_path = os.path.join(params.output_path, "all")
    os.makedirs(all_path)
    for cluster, read_ids in cluster_to_read_ids.iteritems():
        cluster_reads = [reads[read_id] for read_id in read_ids]
        size = len(cluster_reads)
        size_dir = os.path.join(params.output_path, str(size))
        if not os.path.exists(size_dir):
            os.makedirs(size_dir)
        cluster_file_name = cluster + ".fasta"
        cluster_path = os.path.join(size_dir, cluster_file_name)
        with open(cluster_path, "w") as cluster_reads_file:
            SeqIO.write(cluster_reads, cluster_reads_file, "fasta")
        shutil.copy2(cluster_path, os.path.join(all_path, cluster_file_name))
    print "Done writing result"


if __name__ == '__main__':
    main()
