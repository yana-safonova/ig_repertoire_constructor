import sys
from collections import defaultdict

def parse_rcm(filename):
    rcm = {}
    with open(filename) as f:
        for line in f:
            line = [_.strip() for _ in line.strip().split("\t")]
            if len(line) == 2:
                id, cluster = line
                if cluster.strip() == "":
                    cluster = None
            else:
                id, cluster = line[0], None

            rcm[id] = cluster

    return rcm


def dist(s, t):
    if len(s) != len(t):
        print s
        print t
    d = 0
    for i in range(len(s)):
        if s[i] != t[i]:
            d += 1
    return d


def get_barcode(read_name):
    return read_name.split("_UMI:")[-1]


def main():
    # read original barcode groups from amplified_to_orig.rcm
    amplified_to_orig = parse_rcm(sys.argv[1])
    ideal_cluster_to_reads = defaultdict(set)
    for read_name, cluster_id in amplified_to_orig.iteritems():
        # excluding chimeras
        if cluster_id is not None:
            ideal_cluster_to_reads[cluster_id].add(read_name)
    print "Clusters: %d" % len(ideal_cluster_to_reads)

    # for each barcode find all the clusters touching it
    # barcode is found using consensus of all amplicons for each original read separately
    amplified_to_comp = parse_rcm(sys.argv[2])
    barcode_to_clusters = defaultdict(set)
    for cluster_id, read_names in ideal_cluster_to_reads.iteritems():
        non_amplified_to_reads = defaultdict(set)
        for read in read_names:
            non_amplified_to_reads[amplified_to_comp[read]].add(read)
        for original, amplified in non_amplified_to_reads.iteritems():
            barcodes = map(get_barcode, amplified)
            consensus_barcode = list(barcodes[0])
            for i in range(len(consensus_barcode)):
                votes = defaultdict(int)
                for barcode in barcodes:
                    votes[barcode[i]] += 1
                consensus_barcode[i] = max(votes, key = votes.get)
            barcode_to_clusters["".join(consensus_barcode)].add(cluster_id)

    constructed_rcm = parse_rcm(sys.argv[3])
    barcode_to_constructed_clusters = defaultdict(set)
    constructed_cluster_to_reads = defaultdict(set)
    for read, cluster in constructed_rcm.iteritems():
        barcode_to_constructed_clusters[get_barcode(read)].add(cluster)
        constructed_cluster_to_reads[cluster].add(read)

    # for each ideal cluster search for a corresponding constructed cluster
    # count percent of the extra constructed cluster reads and the number of missing ones
    barcodes_to_split = [barcode for barcode in barcode_to_clusters.keys() if len(barcode_to_clusters[barcode]) > 1]
    total_ambiguous_correspondences = 0
    shares_of_extra_reads_in_constructed = []
    shares_of_missing_reads_in_constructed = []
    for barcode in barcodes_to_split:
        for ideal_cluster in barcode_to_clusters[barcode]:
            corresponding_constructed = None
            best_intersection = 0
            ambiguous = False
            for constructed_cluster in barcode_to_constructed_clusters[barcode]:
                reads = constructed_cluster_to_reads[constructed_cluster] & ideal_cluster_to_reads[ideal_cluster]
                if len(reads) > best_intersection:
                    best_intersection = len(reads)
                    corresponding_constructed = constructed_cluster
                    ambiguous = False
                elif len(reads) == best_intersection:
                    ambiguous = True
            if ambiguous:
                total_ambiguous_correspondences += 1

            ideal_reads = set([read for read in ideal_cluster_to_reads[ideal_cluster] if dist(get_barcode(read), barcode) <= 1])
            constructed_reads = set([read for read in constructed_cluster_to_reads[corresponding_constructed] if dist(get_barcode(read), barcode) <= 1])
            shares_of_extra_reads_in_constructed.append((float(len(constructed_reads - ideal_reads)) / len(ideal_reads), len(constructed_reads)))
            shares_of_missing_reads_in_constructed.append((float(len(ideal_reads - constructed_reads)) / len(ideal_reads), len(constructed_reads)))

    shares_of_extra_reads_in_constructed = sorted(shares_of_extra_reads_in_constructed)
    shares_of_missing_reads_in_constructed = sorted(shares_of_missing_reads_in_constructed)

    print "average extra reads in constructed clusters: %f" % (sum(map(lambda x: x[0], shares_of_extra_reads_in_constructed)) / len(shares_of_extra_reads_in_constructed))
    print "median extra reads in constructed clusters: %f" % shares_of_extra_reads_in_constructed[len(shares_of_extra_reads_in_constructed) / 2][0]
    print "average missing reads in constructed clusters: %f" % (sum(map(lambda x: x[0], shares_of_missing_reads_in_constructed)) / len(shares_of_missing_reads_in_constructed))
    print "median missing reads in constructed clusters: %f" % shares_of_missing_reads_in_constructed[len(shares_of_missing_reads_in_constructed) / 2][0]

    with open("missing", "w") as missing, open("extra", "w") as extra:
        for x in shares_of_missing_reads_in_constructed:
            missing.write(str(x[0]) + " " + str(x[1]) + "\n")
        for x in shares_of_extra_reads_in_constructed:
            extra.write(str(x[0]) + " " + str(x[1]) + "\n")

if __name__ == '__main__':
    main()
