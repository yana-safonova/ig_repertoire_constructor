import sys
from Bio import SeqIO
from ig_compress_equal_clusters import parse_cluster_mult
from rcm_utils import read_rcm_list, rcm2rcmint


def merge():
    rcm_file = sys.argv[1]
    id_file = sys.argv[2]
    from_read_file = sys.argv[3]
    to_read_file = sys.argv[4] if sys.argv[4] != "-" else None
    output_rcm_file = sys.argv[5]

    print "reading rcm file from %s" % rcm_file
    rcm = rcm2rcmint(read_rcm_list(rcm_file))

    print "reading id list reported by trie compressor from %s" % id_file
    with open(id_file, "r") as ids:
        id_map = [int(line) for line in ids]
    id_map_set = set(id_map)

    print "reading file with source reads from %s to obtain cluster numbers" % from_read_file
    from_file_cluster_to_pos = dict()
    with open(from_read_file, "r") as reads:
        for i, record in enumerate(SeqIO.parse(reads, "fasta")):
            from_file_cluster_to_pos[int(parse_cluster_mult(record.id)[0])] = i

    to_file_pos_to_cluster = dict()
    if to_read_file is not None:
        print "reading file with destination reads from %s to obtain cluster numbers" % to_read_file
        with open(to_read_file, "r") as reads:
            for i, record in enumerate(SeqIO.parse(reads, "fasta")):
                cluster = int(parse_cluster_mult(record.id)[0])
                to_file_pos_to_cluster[i] = cluster
                id_map_set.remove(i)
        assert not id_map_set
    else:
        current = 0
        for key, value in enumerate(id_map):
            if not value in to_file_pos_to_cluster:
                to_file_pos_to_cluster[value] = current
                current += 1

    result_handle = dict()
    print "merging maps and writing to %s" % output_rcm_file
    with open(output_rcm_file, "w") as output:
        for key, value in rcm:
            handle = to_file_pos_to_cluster[id_map[from_file_cluster_to_pos[value]]] if value is not None else None
            if handle is None:
                output.write("%s\n" % key)
            else:
                output.write("%s\t%d\n" % (key, handle))
            result_handle[key] = handle

    print "checking result"
    # checking if corresponding clusterings are nested
    from collections import defaultdict
    orig_clusters = defaultdict(list)
    for key, value in rcm:
        orig_clusters[value].append(key)
    for orig_handle, ids in orig_clusters.iteritems():
        handle = result_handle[ids[0]]
        for id in ids:
            assert result_handle[id] == handle


if __name__ == "__main__":
    merge()
