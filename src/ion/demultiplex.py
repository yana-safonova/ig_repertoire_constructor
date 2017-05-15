import logging
import sys
import os

home_directory = os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
py_src = os.path.join(home_directory, "py")
sys.path.append(py_src)
print py_src


from ash_python_utils import smart_open


def CreateLogger():
    log = logging.getLogger('find_indels_in_alignment')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    return log


def ParseCommandLineParams(log):
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", type=str, dest="input_file", help="Input file path", required=True)
    parser.add_argument("-m", "--mask", type=str, dest="mask_file", help="Demultiplexing mask file path", required=True)
    parser.add_argument("-o", "--output", type=str, dest="output_file", help="Output file path", required=True)
    parser.add_argument("-b", "--bad", type=str, dest="bad_reads_file", help="Path to a file to fill with reads filtered out", required=True)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    params = parser.parse_args()

    return params


def ReadFastq(reads_path):
    from Bio import SeqIO
    with smart_open(reads_path, "r") as reads_file:
        return list(SeqIO.parse(reads_file, "fastq"))


def WriteFastq(output_path, reads):
    from Bio import SeqIO
    with smart_open(output_path, "w") as reads_file:
        return SeqIO.write(reads, reads_file, "fastq")


def main():
    log = CreateLogger()
    params = ParseCommandLineParams(log)
    reads = ReadFastq(params.input_file)
    # log.info("reading reads from %s" % params.input_file)
    log.info("%d reads read from %s" % (len(reads), params.input_file))

    mask = open(params.mask_file).readline().strip()
    first = mask[:mask.find('N')]
    last = mask[mask.rfind('N') + 1:]
    assert first >= 0 and last >= 0
    print "Searching for %s and %s bounding the barcode." % (first, last)
    # mask = mask.replace('N', '.')
    # barcode_pos = set([p for p in range(len(mask)) if mask[p] == '.'])
    # assert set(mask).issubset(set("ACGT."))

    extracted = 0
    failed = 0
    result = []
    bad = []
    for read in reads:
        # import re
        # m = re.search(mask, str(read.seq))
        # groups = m.groups() if m is not None else []
        # assert len(groups) <= 1
        s = str(read.seq)
        f = s.find(first)
        l = s.find(last)
        # if len(groups) == 0:
        if f < 0 or l < 0 or f + len(first) >= l:
            failed += 1
            bad.append(read)
        else:
            extracted += 1
            barcode = s[f:l + len(last)]
            # barcode = str([ch for i, ch in enumerate(groups[0]) if i in barcode_pos])
            read.id = read.name = read.description = read.id + "|UMI:" + barcode
            result.append(read)

    log.info("Barcode extracted from %d reads. For %d reads extraction failed." % (extracted, failed))

    WriteFastq(params.output_file, result)
    WriteFastq(params.bad_reads_file, bad)
    log.info("Result written to %s. Filtered out reads output to %s" % (params.output_file, params.bad_reads_file))


if __name__ == '__main__':
    main()
