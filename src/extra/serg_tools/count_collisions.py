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


def main():
    rcm = sys.argv[1]
    cluster_to_umi_with_size = dict()
    for line in open(rcm):
        parts = line.split('\t')
        if len(parts) != 2:
            continue
        read_id, cluster_id = parts
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
    print len(collisions)
    print collisions


if __name__ == '__main__':
    main()
