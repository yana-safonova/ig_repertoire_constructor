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
                        help = "path to a file with MIXCR output repertoire (clones.txt)",
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
    records = []
    with open(params.repertoire_path) as input_file:
        header = input_file.readline().split("\t")
        sequence_column = header.index("Clonal sequence(s)")
        size_column = header.index("Clone count")
        id = 0
        for line in input_file:
            if len(line) == 0:
                break
            info = line.split("\t")
            from Bio import Seq
            record = SeqIO.SeqRecord(seq = Seq.Seq(info[sequence_column]), id = "cluster___%d___size___%d" % (id, int(info[size_column])), description = "")
            records.append(record)
            id += 1

    log.info("Read %d reads" % len(records))

    log.info("Writing output")
    with smart_open(params.output_path, "w") as output_file:
        for record in records:
            SeqIO.write(record, output_file, "fasta")


if __name__ == '__main__':
    main()
