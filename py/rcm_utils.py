def read_rcm_map(rcm_file_path):
    rcm = {}
    for line in open(rcm_file_path, "r"):
        key, value = line.strip().split("\t")
        rcm[key.strip()] = value.strip()
    return rcm


def read_rcm_list(rcm_file_path):
    rcm = []
    for line in open(rcm_file_path, "r"):
        key, value = line.strip().split("\t")
        rcm.append((key.strip(), value.strip()))
    return rcm


def rcm2rcmint(rcm):
    # convert values to ints, to be got rid of
    if isinstance(rcm, dict):
        for read in rcm:
            rcm[read] = int(rcm[read])
    else:
        for i, (read, cluster) in enumerate(rcm):
            rcm[i] = (read, int(cluster))
    return rcm


def combine_rcms(rcm1to2_list, ids2_list, rcm2to3_map):
    # Consider two consecutive clusterings of original set of reads with ids id1 to set with ids id2 and then id3
    # The function returns the mapping from original reads to final clusters
    return [(id1, rcm2to3_map[ids2_list[id2]]) for (id1, id2) in rcm1to2_list]
