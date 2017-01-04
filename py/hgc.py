#!/usr/bin/env python2

from test_on_naive_simulation import *
import multiprocessing
from joblib import Parallel, delayed


def run_compressor(input, output):
    run_trie_compressor = os.path.join(igrec_dir, 'build/release/bin/./ig_trie_compressor')
    cmd_line = "%s -i %s -o %s -Toff" % (run_trie_compressor, input, output)
    os.system(cmd_line)


def run_hgc_complexity_estimator(input, output, tau=4):
    run_estimator = os.path.join(igrec_dir, 'build/release/bin/ig_hgc_complexity_estimator')
    cmd_line = "%s -i %s -o %s --tau=%d" % (run_estimator, input, output, tau)
    os.system(cmd_line)


def hgc_estimator(input_file, output_dir, tau=4):
    mkdir_p(output_dir)
    run_compressor(input_file, output_dir + "/compressed.fa.gz")
    run_hgc_complexity_estimator(output_dir + "/compressed.fa.gz",
                                 output_dir + "/compl_stats.txt",
                                 tau=tau)

if __name__ == "__main__":
    datas = ["/Sid/ysafonova/ig_repertoire_constructor/SRR36200%d_igrec/vj_finder/cleaned_reads.fa" % n for n in [74, 75, 98]]
    names = ["SRR74", "SRR75", "SRR98"]

    for data, name in zip(datas, names):
        hgc_estimator(data, name, tau=4)
