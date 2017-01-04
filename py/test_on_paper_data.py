#!/usr/bin/env python2

from test_on_naive_simulation import *
import multiprocessing
from joblib import Parallel, delayed


def test_on_aq_data(name, out_dir, threads=16, do_not_run=False):
    run_and_quast_all(igrec_dir + "/aimquast_test_dataset/%s/input_reads.fa.gz" % name,
                      igrec_dir + "/aimquast_test_dataset/%s/repertoire.fa.gz" % name,
                      igrec_dir + "/aimquast_test_dataset/%s/repertoire.rcm" % name,
                      out_dir=out_dir,
                      threads=threads,
                      do_not_run=do_not_run,
                      do_run_igrec_old=True)


if __name__ == "__main__":
    # test_on_aq_data("test", "test")
    # test_on_aq_data("SIMULATED_1", "SIMULATED_1")
    # test_on_aq_data("SYNTHETIC_1", "SYNTHETIC_1")
    # test_on_aq_data("age3", "age3")
    # test_on_aq_data("flu", "flu_all")

    queue = []
    # queue.append(("test", "test"))
    # queue.append(("SIMULATED_0.5", "SIMULATED_0.5"))
    # queue.append(("SIMULATED_1", "SIMULATED_1"))
    # queue.append(("SIMULATED_2", "SIMULATED_2"))
    queue.append(("SYNTHETIC_0.5", "SYNTHETIC_0.5"))
    queue.append(("SYNTHETIC_1", "SYNTHETIC_1"))
    queue.append(("SYNTHETIC_2", "SYNTHETIC_2"))

    def JOB(a):
        test_on_aq_data(*a)

    n_jobs = 1 if multiprocessing.cpu_count() <= 16 else 3
    Parallel(n_jobs=n_jobs)(delayed(JOB)(q) for q in queue)
