import logging
from umi_tools.count_equal_reads_with_close_umis import smart_open


MAX_READ_DIST = 10
SUM_SIZE_THRESHOLD = 5
MAX_SIZE_THRESHOLD = 3


def CreateLogger():
    log = logging.getLogger('igrec')
    log.setLevel(logging.DEBUG)
    import sys
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    return log


def ParseCommandLineParams(log):
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-r",
                        dest = "reads_path",
                        type = str,
                        help = "path to a file with cleaned original reads",
                        required = True)
    parser.add_argument("-m",
                        dest = "rcm_path",
                        type = str,
                        help = "path to a file with final repertoire rcm",
                        required = True)
    parser.add_argument("-u",
                        dest = "umi_path",
                        type = str,
                        help = "path to a file with extracted barcodes",
                        required = True)
    parser.add_argument("-o",
                        dest = "output_path",
                        type = str,
                        help = "path to an output file",
                        required = True)

    params = parser.parse_args()

    return params


def ReadInput(log, params):
    from Bio import SeqIO
    log.info("Reading reads from %s" % params.reads_path)
    with smart_open(params.reads_path) as reads_file:
        read_id_to_read = dict([(record.id, record.seq) for record in SeqIO.parse(reads_file, "fasta")])
    log.info("Read %d reads" % len(read_id_to_read))

    log.info("Reading rcm from %s" % params.rcm_path)
    from rcm_utils import read_rcm_list
    rcm = read_rcm_list(params.rcm_path)
    log.info("Read %d mappings" % len(rcm))

    log.info("Reading UMIs from %s" % params.umi_path)
    with smart_open(params.umi_path) as umis_file:
        umis = list(SeqIO.parse(umis_file, "fastq"))
    log.info("Read %d UMIs" % len(umis))

    return read_id_to_read, rcm, umis


def ReportCorrectedUmiErrors(log, read_id_to_read, rcm, umis, output_path):
    read_id_to_umi = dict()
    umi_to_read_ids = dict()
    for record in umis:
        read_id_to_umi[record.id] = str(record.seq)
        if record.seq not in umi_to_read_ids:
            umi_to_read_ids[record.seq] = []
        umi_to_read_ids[record.seq].append(record.id)
    assert len(umis) == len(read_id_to_umi)

    cluster_to_read_id = dict()
    for read, cluster in rcm:
        if cluster not in cluster_to_read_id:
            cluster_to_read_id[cluster] = []
        cluster_to_read_id[cluster].append(read)

    output_file = smart_open(output_path, "w")
    adjacent_pairs = 0
    large_pairs = 0
    for cluster, ids in cluster_to_read_id.iteritems():
        # umis = set([str(read_id_to_umi[read_id]) for read_id in ids])
        cluster_umis = set()
        for read_id in ids:
            cluster_umis.add(read_id_to_umi[read_id])
        for umi1 in cluster_umis:
            for umi2 in cluster_umis:
                if umi1 >= umi2:
                    continue
                from umi_tools.count_equal_reads_with_close_umis import dist
                if dist(umi1, umi2) == 1:
                    adjacent_pairs += 1
                    reads1 = set()
                    reads2 = set()
                    for read1 in umi_to_read_ids[umi1]:
                        for read2 in umi_to_read_ids[umi2]:
                            if dist(read1, read2, MAX_READ_DIST + 1) <= MAX_READ_DIST:
                                reads1.add(read1)
                                reads2.add(read2)
                    if len(reads1) + len(reads2) < SUM_SIZE_THRESHOLD or max(len(reads1), len(reads2)) < MAX_SIZE_THRESHOLD:
                        continue
                    large_pairs += 1
                    output_file.write("New pair: %s %s\n" % (umi1, umi2))
                    output_file.write("%d reads from first + %d reads from second\n" % (len(reads1), len(reads2)))
                    for read1 in reads1:
                        output_file.write("%s\n" % read_id_to_read[read1])
                    output_file.write("\n")
                    for read2 in reads2:
                        output_file.write("%s\n" % read_id_to_read[read2])
                    output_file.write("-----------------")

    log.info("Total %d of adjacent pairs found" % adjacent_pairs)
    log.info("Total %d of them have at least %d in total and at least %d in the largest" % (large_pairs, SUM_SIZE_THRESHOLD, MAX_SIZE_THRESHOLD))


def main():
    log = CreateLogger()
    params = ParseCommandLineParams(log)
    read_id_to_read, rcm, umis = ReadInput(log, params)
    ReportCorrectedUmiErrors(log, read_id_to_read, rcm, umis, params.output_path)


if __name__ == '__main__':
    main()
