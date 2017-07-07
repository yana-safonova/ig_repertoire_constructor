import sys
from collections import defaultdict


CLUSTER_SIZE_THRESHOLD = 0
MIN_CLUSTER_SHARE = 5


def dist(s, t):
    if len(s) != len(t):
        print s
        print t
    d = 0
    for i in range(len(s)):
        if s[i] != t[i]:
            d += 1
    return d


def get_vj_hits(alignment_path):
    with open(alignment_path) as alignment:
        columns = alignment.readline().strip().split("\t")
        id_index = columns.index("Read_name")
        v_index = columns.index("V_hit")
        j_index = columns.index("J_hit")
        vj_hits = dict()
        for line in alignment:
            items = line.strip().split("\t")
            vj_hits[items[id_index]] = items[v_index], items[j_index]
    return vj_hits


def main():
    rcm = sys.argv[1]
    cluster_to_umi_with_size = dict()
    read_to_cluster = dict()
    for line in open(rcm):
        parts = line.strip().split('\t')
        if len(parts) != 2:
            continue
        read_id, cluster_id = parts
        read_to_cluster[read_id] = cluster_id
        umi = read_id.split('UMI:')[1].split(':')[0]
        if cluster_id not in cluster_to_umi_with_size:
            cluster_to_umi_with_size[cluster_id] = defaultdict(int)
        cluster_to_umi_with_size[cluster_id][umi] += 1

    for cluster_id in cluster_to_umi_with_size.keys():
        umi_to_size = cluster_to_umi_with_size[cluster_id]
        size = sum(umi_to_size.values())
        if size < CLUSTER_SIZE_THRESHOLD:
            del cluster_to_umi_with_size[cluster_id]
            continue
        gone = set()
        while not set(umi_to_size.keys()).issubset(gone):
            left = set(umi_to_size.keys()) - gone
            min_umi = min(left, key=umi_to_size.get)
            for umi in umi_to_size.keys():
                if dist(min_umi, umi) == 1 and umi_to_size[min_umi] > umi_to_size[umi]:
                    del umi_to_size[umi]
            gone.add(min_umi)

    umi_to_cluster_shares = dict()
    for cluster_id, umi_to_size in cluster_to_umi_with_size.iteritems():
        for umi in umi_to_size:
            if umi not in umi_to_cluster_shares:
                umi_to_cluster_shares[umi] = defaultdict(int)
            umi_to_cluster_shares[umi][cluster_id] = umi_to_size[umi]

    for umi in umi_to_cluster_shares.keys():
        shares = sorted(umi_to_cluster_shares[umi].values())
        if len(shares) < 2 or shares[-2] < MIN_CLUSTER_SHARE:
            del umi_to_cluster_shares[umi]


    collisions = [(umi, shares) for umi, shares in umi_to_cluster_shares.iteritems()]
    print "collisions: %d" % len(collisions)
    print collisions


    vj_hits = get_vj_hits(sys.argv[2])
    cluster_v_hits = defaultdict(set)
    cluster_j_hits = defaultdict(set)
    for read_id, cluster_id in read_to_cluster.iteritems():
        v_hit, j_hit = vj_hits[read_id]
        cluster_v_hits[cluster_id].add(v_hit)
        cluster_j_hits[cluster_id].add(j_hit)


    for umi in umi_to_cluster_shares.keys():
        shares = sorted(umi_to_cluster_shares[umi].values())
        cluster1 = shares[-1]
        cluster2 = shares[-2]
        if not cluster_v_hits[cluster1].isdisjoint(cluster_v_hits[cluster2]):
            continue
        if not cluster_j_hits[cluster1].isdisjoint(cluster_j_hits[cluster2]):
            continue
        del umi_to_cluster_shares[umi]

    strict_collisions = [(umi, shares) for umi, shares in umi_to_cluster_shares.iteritems()]
    print "strict collisions: %d" % len(strict_collisions)
    print strict_collisions


if __name__ == '__main__':
    main()
