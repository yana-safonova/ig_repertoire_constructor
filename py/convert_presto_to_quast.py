import logging
from umi_tools.count_equal_reads_with_close_umis import smart_open


def CreateLogger():
    log = logging.getLogger('convert_presto_to_quast')
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
                        dest = "repertoire_path",
                        type = str,
                        help = "path to a file with pRESTO output repertoire",
                        required = True)
    parser.add_argument("-o",
                        dest = "output_path",
                        type = str,
                        help = "path to an output file",
                        required = True)

    params = parser.parse_args()

    return params


def main():
    log = CreateLogger()
    params = ParseCommandLineParams(log)
    from Bio import SeqIO
    log.info("Reading reads from %s" % params.repertoire_path)
    with smart_open(params.repertoire_path) as reads_file:
        records = list(SeqIO.parse(reads_file, "fasta"))
    log.info("Read %d reads" % len(records))

    log.info("Converting records")
    cnt = 0
    for record in records:
        id = record.id
        import re
        match = re.match(r"^.*CONSCOUNT=(\d+).*$", id)
        record.description = ""
        record.id = "cluster___%d___size___%d" % (cnt, int(match.groups()[0]))
        cnt += 1

    log.info("Writing output")
    with smart_open(params.output_path, "w") as output_file:
        for record in records:
            SeqIO.write(record, output_file, "fasta")


if __name__ == '__main__':
    main()
