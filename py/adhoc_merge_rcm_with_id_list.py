import sys
from rcm_utils import read_rcm_list, rcm2rcmint


def merge():
    rcm_file = sys.argv[1]
    id_file = sys.argv[2]
    output_rcm_file = sys.argv[3]
    print "reading rcm file from %s" % rcm_file
    rcm = rcm2rcmint(read_rcm_list(rcm_file))
    print "reading id list reported by trie compressor from %s" % id_file
    with open(id_file, "r") as ids:
        id_map = [int(line) for line in ids]
    print "merging maps and writing to %s" % output_rcm_file
    with open(output_rcm_file, "w") as output:
        for key, value in rcm:
            output.write("%s\t%d\n" % (key, id_map[value]))


if __name__ == "__main__":
    merge()