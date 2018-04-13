#!/usr/bin/env python2

from test_on_naive_simulation import *
import multiprocessing
from joblib import Parallel, delayed


def test_on_aq_data(name, out_dir, threads=16, do_not_run=False):
    run_and_quast_all(igrec_dir + "/test_dataset/igquast/%s/input_reads.fa.gz" % name,
                      igrec_dir + "/test_dataset/igquast/%s/repertoire.fa.gz" % name,
                      igrec_dir + "/test_dataset/igquast/%s/repertoire.rcm" % name,
                      out_dir=out_dir,
                      threads=threads,
                      do_not_run=do_not_run,
                      # do_run_igrec_old=True)
                      do_run_igrec_old=False)


if __name__ == "__main__":
    queue = []
    queue.append(("test", "test"))
    # queue.append(("age3", "REAL"))
    # queue.append(("flu", "FLU"))
    # queue.append(("SIMULATED_0.5", "SIMULATED_0.5"))
    # queue.append(("SIMULATED_1", "SIMULATED_1"))
    # queue.append(("SIMULATED_2", "SIMULATED_2"))
    queue.append(("donor_1_S1", "donor_1_S1"))
    queue.append(("donor_2_S1", "donor_2_S1"))
    queue.append(("donor_1_S2", "donor_1_S2"))
    queue.append(("donor_2_S2", "donor_2_S2"))
    # queue.append(("donor_1_S1_stage1", "donor_1_S1_stage1"))
    # queue.append(("donor_1_S1_stage3", "donor_1_S1_stage3"))
    # queue.append(("SYNTHETIC_0.5", "SYNTHETIC_0.5"))
    # queue.append(("SYNTHETIC_1", "SYNTHETIC_1"))
    # queue.append(("SYNTHETIC_2", "SYNTHETIC_2"))
    queue.append(("age3_chu", "REAL_CHU"))

    stern = "Stern_SRR1383326 Stern_SRR1383447 Stern_SRR1383448 Stern_SRR1383449 Stern_SRR1383450 Stern_SRR1383451 Stern_SRR1383452 Stern_SRR1383453 Stern_SRR1383454 Stern_SRR1383455 Stern_SRR1383456 Stern_SRR1383457 Stern_SRR1383458 Stern_SRR1383459 \
        Stern_SRR1383460 Stern_SRR1383461 Stern_SRR1383462 Stern_SRR1383463 Stern_SRR1383464 Stern_SRR1383465 Stern_SRR1383466 Stern_SRR1383467 Stern_SRR1383468 Stern_SRR1383469 Stern_SRR1383470 Stern_SRR1383471 Stern_SRR1383472 \
        Stern_SRR1383473 Stern_SRR1383474 Stern_SRR1383475 Stern_SRR1383476 Stern_SRR1383477".split()
    # stern = "Stern_SRR1383326 Stern_SRR1383447 Stern_SRR1383448 Stern_SRR1383449 Stern_SRR1383450 Stern_SRR1383451 Stern_SRR1383452 Stern_SRR1383453 Stern_SRR1383454 Stern_SRR1383455 Stern_SRR1383456 Stern_SRR1383457 Stern_SRR1383458 Stern_SRR1383459 \
    #     Stern_SRR1383460".split()
    for s in stern:
        queue.append((s, s))

    datasets = ["naive_no_pcr", "naive_with_pcr", "plasma_no_pcr", "plasma_with_pcr"]
    for dataset in datasets:
        queue.append((dataset, dataset))

    for i in range(1, 10):
        dataset = "AGE%d" % i
        queue.append((dataset, dataset))

    def JOB(a):
        test_on_aq_data(*a)

    n_jobs = 1 if multiprocessing.cpu_count() <= 16 else 3
    Parallel(n_jobs=n_jobs)(delayed(JOB)(q) for q in queue)
    # for q in queue:
    #     JOB(q)
