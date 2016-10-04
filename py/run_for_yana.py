#!/usr/bin/env python2

from simulate import *

if __name__ == "__main__":
    data1 = "/Sid/ysafonova/ig_repertoire_constructor/BSSE_naive/vj_finder/cleaned_reads.fa"
    data2 = "/Sid/ysafonova/ig_repertoire_constructor/BSSE_asc/vj_finder/cleaned_reads.fa"
    data3 = "/Sid/ysafonova/ig_repertoire_constructor/BSSE_plasma/vj_finder/cleaned_reads.fa"

    datas = [data1, data2, data3]
    names = ["naive", "asc", "plasma"]

    for data, name in zip(datas, names):
        run_mixcr(data, "mice_out_mixcr_" + name, species="mmu", remove_tmp=False)
        run_presto(data, "mice_out_presto_" + name)
