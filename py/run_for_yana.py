#!/usr/bin/env python2

from simulate import *

import sys

if __name__ == "__main__":
    # run_mixcr2("test_dataset/test7.fastq", "test7_mixcr2", remove_tmp=False)
    # run_mixcr("test_dataset/test7.fastq", "test7")

    # sys.exit()

    # data1 = "/Sid/ysafonova/ig_repertoire_constructor/BSSE_naive/vj_finder/cleaned_reads.fa"
    # data2 = "/Sid/ysafonova/ig_repertoire_constructor/BSSE_asc/vj_finder/cleaned_reads.fa"
    # data3 = "/Sid/ysafonova/ig_repertoire_constructor/BSSE_plasma/vj_finder/cleaned_reads.fa"

    # datas = [data1, data2, data3]
    # names = ["naive", "asc", "plasma"]
    #
    # data1 = "/Sid/ysafonova/ig_repertoire_constructor/SRR3620074.fastq"
    # data2 = "/Sid/ysafonova/ig_repertoire_constructor/SRR3620075.fastq"
    # data3 = "/Sid/ysafonova/ig_repertoire_constructor/SRR3620098.fastq"

    datas = ["/Sid/ysafonova/ig_repertoire_constructor/SRR36200%d_igrec/vj_finder/cleaned_reads.fa" % n for n in [74, 75, 98]]
    names = ["SRR74", "SRR75", "SRR98"]

    for data, name in zip(datas, names):
        run_mixcr2(data, "mixcr_CDR1_CDR3_" + name, species="hsa", remove_tmp=True,
                   region_from="FR1End", region_to="FR4Begin")
        run_mixcr2(data, "mixcr_FR1_CDR3_" + name, species="hsa", remove_tmp=True,
                   region_from="FR1Begin", region_to="FR4Begin")
        run_mixcr2(data, "mixcr_FR1_FR4_" + name, species="hsa", remove_tmp=True,
                   region_from="FR1Begin", region_to="FR4End")
        run_presto(data, "presto_" + name)
        run_igrec(data, "igrec_" + name, additional_args=" --no-alignment")
        # data_after_vjf = "igrec_" + name + "/vj_finder/cleaned_reads.fa"
        # run_mixcr2(data_after_vjf, "mixcr_after_vjf_CDR1_CDR3_" + name, species="hsa", remove_tmp=True,
        #            region_from="FR1End", region_to="FR4Begin")
        # run_mixcr2(data_after_vjf, "mixcr_after_vjf_FR1_CDR3_" + name, species="hsa", remove_tmp=True,
        #            region_from="FR1Begin", region_to="FR4Begin")
        # run_mixcr2(data_after_vjf, "mixcr_after_vjf_FR1_FR4_" + name, species="hsa", remove_tmp=True,
        #            region_from="FR1Begin", region_to="FR4End")
