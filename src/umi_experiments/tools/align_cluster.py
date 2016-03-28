import logging
import os
import sys
from collections import defaultdict
from Bio import SeqIO


def CreateLogger():
    log = logging.getLogger('align_cluster')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    return log


def ParseCommandLineParams(log):
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-f", "--final-repertoire-rcm", type=str,   dest="final_rcm_path",          help="Path to final repertoire rcm",            required=True)
    parser.add_argument("-i", "--inter-repertoire-rcm", type=str,   dest="inter_rcm_path",          help="Path to intermediate repertoire rcm",     required=True)
    parser.add_argument("-c", "--inter-clusters",       type=str,   dest="inter_fastq_path",        help="Path to intermediate repertoire fastq",   required=True)
    parser.add_argument("-r", "--reads",                type=str,   dest="reads_path",              help="Path to fastq file with reads",           required=True)
    parser.add_argument("-o", "--output",               type=str,   dest="output_dir",              help="Path to output directory",                required=True)
    # parser.add_argument("-c", "--cluster-id",           type=int,   dest="cluster_id",  default=0,  help="Id of the cluster to analyze",        required=False)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    params = parser.parse_args()

    return params


def FindAbundantCluster(log, rcm_path):
    with open(rcm_path, "r") as rcm:
        for line in rcm:
            id_and_cluster = line[:-1].split("\t")
            size = int(id_and_cluster[0].split("___")[-1])
            if size >= 5:
                log.info("Found cluster '%s', %s of size %d", id_and_cluster[0], id_and_cluster[1], size)
                return int(id_and_cluster[1])
    assert False


# Relies that rcm contains reads in order
def ExtractReadsFromCluster(log, rcm_path, cluster_id, read_id_to_idx):
    log.info("Extracting cluster reads from %s", rcm_path)
    read_ids = []
    with open(rcm_path, "r") as rcm:
        for line in rcm:
            id_and_cluster = line.split("\t")
            if int(id_and_cluster[1]) == cluster_id:
                log.debug("Read id: '%s', its cluster: '%d', cluster id: %d", id_and_cluster[0], int(id_and_cluster[1]), cluster_id)
                read_ids.append((read_id_to_idx[id_and_cluster[0]], id_and_cluster[0]))
    log.info("Extracted %d reads from cluster %d (rcm %s)", len(read_ids), cluster_id, rcm_path)
    return read_ids


def ReadRecords(log, reads_path):
    with open(reads_path, "r") as reads_file:
        records_list = list(SeqIO.parse(reads_file, "fastq"))
        id_to_record = SeqIO.to_dict(records_list)
        id_to_idx = defaultdict()
        for i in range(len(records_list)):
            id_to_idx[records_list[i].id] = i
        log.info("Read %d (%d) records from %s", len(id_to_record), len(records_list), reads_path)
        return id_to_idx, id_to_record


def GenFastqNClustal(log, reads, dir, name):
    path = os.path.join(dir, name)
    log.info("Writing subcluster of %d reads to %s", len(reads), path)
    with open(path) as cluster_file:
        SeqIO.write(reads, cluster_file, "fastq")
    clustal_cmd = "clustalw -infile=%s" % path
    log.info("Running '%s'", clustal_cmd)
    os.system(clustal_cmd)


def main():
    log = CreateLogger()
    params = ParseCommandLineParams(log)
    final_cluster_id = FindAbundantCluster(log, params.final_rcm_path)
    inter_id_to_idx, inter_id_to_record = ReadRecords(log, params.inter_fastq_path)
    clusters = ExtractReadsFromCluster(log, params.final_rcm_path, final_cluster_id, inter_id_to_idx)
    read_to_idx, read_id_to_record = ReadRecords(log, params.reads_path)
    all_reads = []
    for cluster in clusters:
        cluster_read_info = ExtractReadsFromCluster(log, params.inter_rcm_path, cluster[0], read_to_idx)
        cluster_reads = []
        for read_info in cluster_read_info:
            cluster_reads.append(read_id_to_record[read_info[1]])
        all_reads.extend(cluster_reads)
        GenFastqNClustal(log, cluster_reads, params.output_dir, str(cluster[0]) + ".fastq")

    GenFastqNClustal(log, all_reads, params.output_dir, "all.fastq")


if __name__ == '__main__':
    main()
