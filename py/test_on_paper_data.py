#!/usr/bin/env python2

from test_on_naive_simulation import *


def test_on_aq_data(name, out_dir, threads=16, do_not_run=False):
    run_and_quast_all(igrec_dir + "/aimquast_test_dataset/%s/input_reads.fa.gz" % name,
                      igrec_dir + "/aimquast_test_dataset/%s/repertoire.fa.gz" % name,
                      igrec_dir + "/aimquast_test_dataset/%s/repertoire.rcm" % name,
                      out_dir=out_dir,
                      threads=threads,
                      do_not_run=do_not_run,
                      do_run_igrec_old=True)


if __name__ == "__main__":
    test_on_aq_data("age3", "age3_all")
    # test_on_aq_data("flu", "flu_all")
