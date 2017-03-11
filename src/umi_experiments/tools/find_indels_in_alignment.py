import logging
import os
import sys
from collections import defaultdict
from Bio import AlignIO

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
    parser.add_argument("-i", "--input", type=str, dest="input_dir", help="Input directory path", required=True)
    parser.add_argument("-o", "--output", type=str, dest="output_file", help="Output statistics file path", required=True)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    params = parser.parse_args()

    return params


def HammingDist(first, second):
    r = 0
    for c, d in zip(first, second):
        if c != d:
            r += 1
    return r

def GetReadsInLargestCluster(alignment, num_columns):
    center = list(alignment[0])
    for i in range(num_columns):
        freq = defaultdict(int)
        for read in alignment:
            freq[read[i]] += 1
        max_freq = max(freq.values())
        # for k, v in freq.iteritems():
        #     if v == max_freq:
        #         center[i] = k
        #         break
        center[i] = max(freq, key = freq.get)
    reads = []
    for read in alignment:
        if HammingDist(read, center) <= 10:
            reads.append(read)
    return reads

def main():
    log = CreateLogger()
    params = ParseCommandLineParams(log)

    indels = defaultdict(int)
    indel_lengths = dict()
    max_columns = 0

    files_num = len(os.listdir(params.input_dir))
    done = 0
    for file in os.listdir(params.input_dir):
        alignment = AlignIO.read(file, "clustal")
        num_columns = len(alignment[0])
        max_columns = max(max_columns, num_columns)
        reads = GetReadsInLargestCluster(alignment, num_columns)
        printed = False
        for read in reads:
            l = 0
            while l < num_columns:
                while l < num_columns and read[l] != '-':
                    l += 1
                if l == num_columns:
                    continue
                if l == 0 and not printed:
                    print file
                    print 'selected:', len(reads)
                    for tmp in reads:
                        print str(tmp.seq)
                    print 'all:', len(alignment)
                    for tmp in alignment:
                        print str(tmp.seq)
                    print
                    printed = True
                r = l + 1
                while r < num_columns and read[r] == '-':
                    r += 1
                if r >= num_columns:
                    break
                for j in range(l, r):
                    indels[j] += 1
                    if not j in indel_lengths:
                        indel_lengths[j] = defaultdict(int)
                    indel_lengths[j][r - l] += 1
                l = r

        done += 1
        print done, "of", files_num

    out = open(params.output_file, "w")
    out.write("pos\tindels\tindel len\tcount\tindel len\tcount\n")
    for j in range(max_columns):
        out.write(`j` + "\t" + `indels[j]`)
        if j in indel_lengths:
            for l in sorted(indel_lengths[j]):
                out.write("\t" + `l` + "\t" + `indel_lengths[j][l]`)
        out.write("\n")
    out.close()

if __name__ == '__main__':
    main()
