def read_rcm_map(rcm_file_path):
    rcm = {}
    for line in open(rcm_file_path, "r"):
        key_value = line.strip().split("\t")
        if len(key_value) == 2:
            key, value = key_value
        elif len(key_value) == 1:
            key = key_value[0]
            value = None
        else:
            raise "Bad RCM line '%s'" % line
        rcm[key.strip()] = value.strip() if value is not None else value
    return rcm


def read_rcm_list(rcm_file_path):
    rcm = []
    for line in open(rcm_file_path, "r"):
        key_value = line.strip().split("\t")
        if len(key_value) == 2:
            key, value = key_value
        elif len(key_value) == 1:
            key = key_value[0]
            value = None
        else:
            raise "Bad RCM line '%s'" % line
        rcm.append((key.strip(), value.strip() if value is not None else value))
    return rcm


def rcm2rcmint(rcm):
    # convert values to ints, to be got rid of
    if isinstance(rcm, dict):
        for read in rcm:
            rcm[read] = int(rcm[read]) if rcm[read] is not None else None
    else:
        for i, (read, cluster) in enumerate(rcm):
            rcm[i] = (read, int(cluster) if cluster is not None else None)
    return rcm


def combine_rcms(rcm1to2_list, ids2_list, rcm2to3_map):
    # Consider two consecutive clusterings of original set of reads with ids id1 to set with ids id2 and then id3
    # The function returns the mapping from original reads to final clusters
    return [(id1, rcm2to3_map[ids2_list[id2]]) for (id1, id2) in rcm1to2_list]
