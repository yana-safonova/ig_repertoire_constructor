import os
from rcm_utils import read_rcm_list, rcm2rcmint


def merge():
    print os.getcwd()
    print "reading intermed rcm"
    inter_rcm = rcm2rcmint(read_rcm_list("umi_clustering/intermediate_repertoire.rcm"))
    print "reading final rcm"
    igrec_rcm_list = read_rcm_list("final_repertoire/final_repertoire.rcm")
    cluster_id_list = {}
    print "mapping"
    for mapping in igrec_rcm_list:
        import re
        m = re.match(r"^intermediate_cluster___(\d+)___size___(\d+)$", mapping[0])
        cluster_id_list[int(m.groups()[0])] = mapping[0]
    igrec_rcm_map = dict(igrec_rcm_list)
    print "writing"
    with open("final_repertoire/final_repertoire_merged.rcm", "w") as fout:
        for original_cleaned, cluster_idx in inter_rcm:
            fout.write(original_cleaned + "\t" + (igrec_rcm_map[cluster_id_list[cluster_idx]] if (cluster_idx in cluster_id_list) else "") + "\n")
    print "done"


if __name__ == "__main__":
    merge()