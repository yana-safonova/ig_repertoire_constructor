import logging


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
                        help = "path to a file with cleaned input reads",
                        required = True)
    parser.add_argument("-u",
                        dest = "umi_path",
                        type = str,
                        help = "path to a file with extracted barcodes",
                        required = True)

    params = parser.parse_args()

    return params


def smart_open(filename, mode="r"):
    import gzip
    import re

    if "w" in mode:
        MODE = "w"
    elif "a" in mode:
        MODE = "a"
    else:
        MODE = "r"

    if re.match(r"^.*\.gz$", filename):
        assert(MODE != "a")
        fh = gzip.open(filename, mode=MODE)
    else:
        fh = open(filename, mode=mode)
    return fh


def ReadInput(log, params):
    from Bio import SeqIO
    log.info("Reading reads")
    with smart_open(params.reads_path) as reads_file:
        reads = [record.seq for record in SeqIO.parse(reads_file, "fasta")]
    log.info("Reading UMIs")
    with smart_open(params.umi_path) as umis_file:
        umis = [record.seq for record in SeqIO.parse(umis_file, "fasta")]
    return reads, umis


def dist(umi1, umi2, max = float("inf")):
    result = 0
    for i in range(min(len(umi1), len(umi2))):
        result += 1 if umi1[i] != umi2[i] else 0
        if result > max:
            return result
    return result# + abs(len(umi1) - len(umi2))


def min_pair(umi1, umi2):
    return (umi1, umi2) if umi1 < umi2 else (umi2, umi1)


def CountEqual(log, reads, umis):
    log.info("Counting stats")
    read_to_umis = dict()
    for read, umi in zip(reads, umis):
        if read not in read_to_umis:
            read_to_umis[read] = set()
        read_to_umis[read].add(umi)

    result = set()
    for read, umi_set in read_to_umis.iteritems():
        for umi1 in umi_set:
            for umi2 in umi_set:
                if dist(umi1, umi2) == 1:
                    result.add(min_pair(umi1, umi2))

    log.info("%d pairs of adjacent UMIs with identical reads found" % len(result))


def main():
    log = CreateLogger()
    params = ParseCommandLineParams(log)
    reads, umis = ReadInput(log, params)
    CountEqual(log, reads, umis)


if __name__ == '__main__':
    main()
