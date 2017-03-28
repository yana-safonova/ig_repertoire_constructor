import os
import sys
from collections import defaultdict

from Bio import SeqIO
current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = os.path.join(current_dir, os.pardir, os.pardir, os.pardir, os.pardir)
sys.path.append(igrec_dir)
sys.path.append(igrec_dir + "/py/")
sys.path.append(igrec_dir + "/py/")
from ig_compress_equal_clusters import parse_cluster_mult
from igquast_impl import parse_rcm
from ash_python_utils import smart_open


def main():
    # read original barcode groups from amplified_to_orig.rcm
    amplified_to_orig = parse_rcm(sys.argv[1])
    id_to_cluster = defaultdict(set)
    for read_name, cluster_id in amplified_to_orig.iteritems():
        # excluding chimeras
        if cluster_id is not None:
            id_to_cluster[cluster_id].add(read_name)
    print "Clusters: %d" % len(id_to_cluster)

    # group them by barcodes (calculated as consensus)
    barcode_to_clusters = defaultdict(set)
    for cluster_id, read_names in id_to_cluster.iteritems():
        barcodes = []
        for read_name in read_names:
            barcodes.append(read_name.split("_UMI:")[1])
        consensus_barcode = list(barcodes[0])
        for i in range(len(consensus_barcode)):
            votes = defaultdict(int)
            for barcode in barcodes:
                votes[barcode[i]] += 1
            consensus_barcode[i] = max(votes, key = votes.get)
        barcode_to_clusters["".join(consensus_barcode)].add(cluster_id)
    print "Unique barcodes: %d" % len(barcode_to_clusters)

    # read final repertoire rcm
    constructed_rcm = parse_rcm(sys.argv[2])

    # check if the reads from groups form separate clusters (count bad ones)
    uncertain = 0
    not_split = 0
    for barcode, clusters in barcode_to_clusters.iteritems():
        if len(clusters) < 2:
            continue
        constructed_clusters = set()
        for cluster in clusters:
            cc_votes = defaultdict(int)
            for read_name in id_to_cluster[cluster]:
                cc_votes[constructed_rcm[read_name]] += 1
            best = max(cc_votes.values())
            largest = [cc for cc in cc_votes.keys() if cc_votes[cc] == best]
            if len(largest) > 1:
                uncertain += 1
            constructed_clusters.update(largest)
        if len(constructed_clusters) < len(clusters):
            not_split += len(clusters) - len(constructed_clusters)

    print "Total %d barcodes are not split" % not_split
    print "%d barcodes have uncertain split. If this number is large, the previous metric is unreliable." % uncertain


if __name__ == '__main__':
    main()
