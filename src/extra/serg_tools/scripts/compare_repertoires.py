import os
import sys
from Bio import SeqIO
current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = os.path.join(current_dir, os.pardir, os.pardir, os.pardir, os.pardir)
sys.path.append(igrec_dir)
sys.path.append(igrec_dir + "/py/")
from ig_compress_equal_clusters import parse_cluster_mult


def read_repertoire(path):
    with open(path) as rep_file:
        records = [record for record in SeqIO.parse(path, "fasta")]
    print len(records), 'records read'
    cluster_to_index = {}
    current = 0
    # wha = 0
    for record in records:
        id, mult = parse_cluster_mult(record.id)
        cluster_to_index[id] = current
        current += 1
    return records, cluster_to_index


def main():
    rep1_path = sys.argv[1]
    rep2_path = sys.argv[2]
    rep1, map1 = read_repertoire(rep1_path)
    rep2, map2 = read_repertoire(rep2_path)
    bad = 0
    only_first = 0
    only_second = 0
    for key, value in map1.iteritems():
        if key not in map2:
            only_first += 1
            print "only first"
            print rep1[value]
        elif rep1[value].seq != rep2[map2[key]].seq:
            bad += 1
            print "not equal:"
            print rep1[value]
            print rep2[map2[key]]
            print len(rep1[value].seq), rep1[value].seq
            print len(rep2[map2[key]].seq), rep2[map2[key]].seq
    for key, value in map2.iteritems():
        if key not in map1:
            only_second += 1
            print "only second"
            print rep2[value]
    # only_second = len(map2) - (len(map1) - only_first)
    print bad, "not equal"
    print only_first, "only first"
    print only_second, "only second"


if __name__ == '__main__':
    main()
