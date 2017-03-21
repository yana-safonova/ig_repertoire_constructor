import logging
from umi_tools.count_equal_reads_with_close_umis import smart_open


MAX_READ_DIST = 5
SUM_SIZE_THRESHOLD = 5
MAX_SIZE_THRESHOLD = 3
LARGE_CLUSTER_SIZE = 5
SIGNIFICANT_CLUSTER_SIZE = 3

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
        sequence = str(record.seq)
        read_id_to_umi[record.id] = sequence
        if sequence not in umi_to_read_ids:
            umi_to_read_ids[sequence] = []
        umi_to_read_ids[sequence].append(record.id)
    assert len(umis) == len(read_id_to_umi)

    cluster_to_read_id = dict()
    for read, cluster in rcm:
        if cluster not in cluster_to_read_id:
            cluster_to_read_id[cluster] = []
        cluster_to_read_id[cluster].append(read)

    output_file = smart_open(output_path, "w")
    adjacent_pairs = 0
    large_pairs = 0
    become_large = 0
    both_significant = 0
    corrected_reads = set()
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
                    read_ids1 = set()
                    read_ids2 = set()
                    for read_id1 in umi_to_read_ids[umi1]:
                        for read_id2 in umi_to_read_ids[umi2]:
                            if dist(read_id_to_read[read_id1], read_id_to_read[read_id2], MAX_READ_DIST) <= MAX_READ_DIST:
                                read_ids1.add(read_id1)
                                read_ids2.add(read_id2)
                    # for _ in range(len(umi_to_read_ids[umi1])):
                    #     for read_id11 in umi_to_read_ids[umi1]:
                    #         if read_id11 not in read_ids1:
                    #             continue
                    #         for read_id12 in umi_to_read_ids[umi1]:
                    #             if read_id12 in read_ids1:
                    #                 continue
                    #             if dist(read_id_to_read[read_id11], read_id_to_read[read_id12], MAX_READ_DIST) <= MAX_READ_DIST:
                    #                 read_ids1.add(read_id12)
                    # for _ in range(len(umi_to_read_ids[umi2])):
                    #     for read_id21 in umi_to_read_ids[umi2]:
                    #         if read_id21 not in read_ids2:
                    #             continue
                    #         for read_id22 in umi_to_read_ids[umi2]:
                    #             if read_id22 in read_ids2:
                    #                 continue
                    #             if dist(read_id_to_read[read_id21], read_id_to_read[read_id22], MAX_READ_DIST) <= MAX_READ_DIST:
                    #                 read_ids2.add(read_id22)

                    if max(len(read_ids1), len(read_ids2)) < LARGE_CLUSTER_SIZE and len(read_ids1) + len(read_ids2) >= LARGE_CLUSTER_SIZE:
                        become_large += 1
                    if len(read_ids1) >= SIGNIFICANT_CLUSTER_SIZE and len(read_ids2) >= SIGNIFICANT_CLUSTER_SIZE:
                        both_significant += 1
                    if len(read_ids1) + len(read_ids2) < SUM_SIZE_THRESHOLD or max(len(read_ids1), len(read_ids2)) < MAX_SIZE_THRESHOLD:
                        continue
                    large_pairs += 1
                    output_file.write("New pair: %s %s\n" % (umi1, umi2))
                    output_file.write("%d reads from first + %d reads from second\n" % (len(read_ids1), len(read_ids2)))
                    for read_id1 in read_ids1:
                        output_file.write(">%s\n%s\n" % (read_id1, read_id_to_read[read_id1]))
                    output_file.write("\n")
                    for read_id2 in read_ids2:
                        output_file.write(">%s\n%s\n" % (read_id2, read_id_to_read[read_id2]))
                    output_file.write("-----------------")

                    not_in_cluster = 0
                    for read_id in read_ids1:
                        if read_id not in ids:
                            not_in_cluster += 1
                    for read_id in read_ids2:
                        if read_id not in ids:
                            not_in_cluster += 1
                    print not_in_cluster, " out of ", len(read_ids1) + len(read_ids2), " are not in the cluster actually!!! Sizes: ", len(read_ids1), len(read_ids2)

                    read_ids = read_ids1 if len(read_ids1) < len(read_ids2) else read_ids2
                    for read_id in read_ids:
                        log.info(">%s\n%s" % (read_id, read_id_to_read[read_id]))
                        # if read_id in corrected_reads:
                        #     print "too bad"
                    corrected_reads.update(read_ids)

    # for read_id in corrected_reads:
    #     log.info(">%s\n%s" % (read_id, read_id_to_read[read_id]))

    log.info("Total %d of adjacent pairs found" % adjacent_pairs)
    log.info("Total %d of them have at least %d in total and at least %d in the largest" % (large_pairs, SUM_SIZE_THRESHOLD, MAX_SIZE_THRESHOLD))
    log.info("%d become large >= %d" % (become_large, LARGE_CLUSTER_SIZE))
    log.info("%d have both at least %d" % (both_significant, SIGNIFICANT_CLUSTER_SIZE))
    log.info("%d total reads corrected" % len(corrected_reads))


def main():
    log = CreateLogger()
    params = ParseCommandLineParams(log)
    read_id_to_read, rcm, umis = ReadInput(log, params)
    ReportCorrectedUmiErrors(log, read_id_to_read, rcm, umis, params.output_path)


if __name__ == '__main__':
    main()
