#!/usr/bin/env python2

from test_on_naive_simulation import *


def run_and_quast_all_three(input_dir, threads=16, do_not_run=False):
    # for i in [1, 2, 3]:
    # for i in [1, 3]:
    for i in [3]:
        run_and_quast_all(input_dir + "/input%d.fa.gz" % i,
                          input_dir + "/repertoire%d.fa.gz" % i,
                          input_dir + "/repertoire%d.rcm" % i,
                          out_dir=input_dir + "/filtering%d" % i,
                          threads=threads,
                          do_not_run=do_not_run)
        # for lam in ["0.5", "1"]:
        #     run_and_quast_all(input_dir + "/input%d_jit%s.fa.gz" % (i, lam),
        #                       input_dir + "/repertoire%d.fa.gz" % i,
        #                       input_dir + "/repertoire%d.rcm" % i,
        #                       out_dir=input_dir + "/filtering%d_jit%s" % (i, lam),
        #                       threads=threads)


if __name__ == "__main__":
    import os.path

    datasets = ["AGE3", "MG91M", "HD09M", "AGE7"]
    # datasets = ["AGE3", "MG91M", "HD09M"]
    datasets = ["AGE3"]
    datasets = []
    datasets = ["FLU_FV_27"]

    datasets = ["AGE3", "AGE7", "MG91M_IGH", "MG91M_IGL", "MG91M_IGK", "FLU_FV_21_IGH", "FLU_FV_21_IGL", "FLU_FL_21_IGK", "FLU_FV_27_IGH", "FLU_FV_27_IGK", "FLU_FV_27_IGL"]
    datasets = ["AGE7", "MG91M_IGH", "MG91M_IGL", "MG91M_IGK", "FLU_FV_21_IGH", "FLU_FV_21_IGL", "FLU_FL_21_IGK", "FLU_FV_27_IGH", "FLU_FV_27_IGK", "FLU_FV_27_IGL"]
    datasets = ["FLU_FV_21_IGH", "FLU_FV_21_IGL", "FLU_FL_21_IGK", "FLU_FV_27_IGH", "FLU_FV_27_IGK", "FLU_FV_27_IGL"]
    datasets = ["MG91M_IGH", "MG91M_IGL", "MG91M_IGK", "HD09M_IGH", "HD09M_IGL", "HD09M_IGK"]
    # MG91M is BAAAAAD!
    datasets = ["AGE3", "AGE7", "FLU_FV_21_IGH", "FLU_FV_21_IGL", "FLU_FL_21_IGK", "FLU_FV_27_IGH", "FLU_FV_27_IGK", "FLU_FV_27_IGL"]
    datasets = ["AGE3", "FLU_FV_21_IGH", "FLU_FV_21_IGL", "FLU_FL_21_IGK", "FLU_FV_27_IGH", "FLU_FV_27_IGK", "FLU_FV_27_IGL", "AGE7"]
    datasets = ["AGE3", "FLU_FV_21_IGH", "FLU_FV_21_IGL", "FLU_FL_21_IGK", "FLU_FV_27_IGH", "FLU_FV_27_IGK", "FLU_FV_27_IGL", "AGE7"]

    datasets = ["AGE1", "AGE2", "AGE3", "AGE4", "AGE5", "AGE6", "AGE7", "AGE8", "AGE9"] + \
               ["FLU_FV_21_IGH", "FLU_FV_21_IGL", "FLU_FL_21_IGK", "FLU_FV_27_IGH", "FLU_FV_27_IGK", "FLU_FV_27_IGL"] + \
               ["HD09M_IGH", "HD09M_IGL", "HD09M_IGK"]

    datasets = ["FLU_FV_21_IGH", "FLU_FV_21_IGL", "FLU_FL_21_IGK"]
    datasets = ["AGE3", "FLU_FV_21_IGH", "AGE7"]

    def JOB(dataset):
        ddir = igrec_dir + "/src/extra/ig_quast_tool/" + dataset
        if os.path.isfile(ddir + "/input1.fa.gz"):
            run_and_quast_all_three(ddir, do_not_run=False)
    import multiprocessing
    n_jobs = 1 if multiprocessing.cpu_count() <= 16 else 3
    Parallel(n_jobs=n_jobs)(delayed(JOB)(dataset) for dataset in datasets)
