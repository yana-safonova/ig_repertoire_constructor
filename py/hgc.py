#!/usr/bin/env python2

from test_on_naive_simulation import *
import multiprocessing
from joblib import Parallel, delayed


def run_compressor(input, output):
    run_trie_compressor = os.path.join(igrec_dir, 'build/release/bin/./ig_trie_compressor')
    cmd_line = "%s -i %s -o %s -Toff" % (run_trie_compressor, input, output)
    os.system(cmd_line)


def run_hgc_complexity_estimator(input, output, tau=4, kstep=5):
    run_estimator = os.path.join(igrec_dir, 'build/release/bin/ig_hgc_complexity_estimator')
    cmd_line = "%s -i %s -o %s --tau=%d --k-step=%d" % (run_estimator, input, output, tau, kstep)
    os.system(cmd_line)


def hgc_estimator(input_file, output_dir, tau=4):
    mkdir_p(output_dir)
    run_compressor(input_file, output_dir + "/compressed.fa.gz")
    run_hgc_complexity_estimator(output_dir + "/compressed.fa.gz",
                                 output_dir + "/compl_stats.txt",
                                 tau=tau)


def hgc_estimator_all(input_file, output_dir):
    mkdir_p(output_dir)
    run_compressor(input_file, output_dir + "/compressed.fa.gz")
    for tau in [1, 2, 3, 4]:
        run_hgc_complexity_estimator(output_dir + "/compressed.fa.gz",
                                     output_dir + "/compl_tau%d.txt" % tau,
                                     tau=tau,
                                     kstep=1)

if __name__ == "__main__":
    datas = []
    names = []
    # datas += ["/Sid/ysafonova/ig_repertoire_constructor/SRR36200%d_igrec/vj_finder/cleaned_reads.fa" % n for n in [74, 75, 98]]
    # names += ["SRR74", "SRR75", "SRR98"]

    # datas += ["/home/ashlemov/Git/ig_repertoire_constructor/SRA_performance_test/cleaned_reads/SRR138346%d.fa.gz" % n for n in xrange(6)]
    # names += ["SRR0%d" % n for n in xrange(6)]

    # datas += ["/home/ashlemov/Git/ig_repertoire_constructor/src/extra/ig_quast_tool/AGE%d/input1.fa.gz" % n for n in xrange(1, 10)]
    # names += ["AGE_%d" % n for n in xrange(1, 10)]

    # datas += ["/home/ashlemov/Git/ig_repertoire_constructor/src/extra/ig_quast_tool/AGE%d/input1.fa.gz" % n for n in [9]]
    # names += ["AGE_%d" % n for n in [9]]
    #
    # datas += ["/Nancy/data/input/ImmunoSeq/roche_datasets/%d_SAM133069%d/merged_reads/%d_SAM133069%d.cleaned.fastq" % (n + 1, n + 69, n + 1, n + 69) for n in xrange(8)]
    # names += ["ROCHE_%d" % (n + 1) for n in xrange(8)]
    # datas += ["/home/ashlemov/Git/ig_repertoire_constructor/test_dataset/igquast/age3/input_reads.fa.gz"]
    # names += ["REAL"]
    # datas += ["/home/ashlemov/Git/ig_repertoire_constructor/test_dataset/igquast/SIMULATED_1/input_reads.fa.gz"]
    # names += ["SIMULATED_1"]
    # datas += ["/home/ashlemov/Git/ig_repertoire_constructor/test_dataset/igquast/SYNTHETIC_1/input_reads.fa.gz"]
    # names += ["SYNTHETIC_1"]
    datas += ["/home/ashlemov/Git/ig_repertoire_constructor/test_dataset/igquast/SIMULATED_0.5/input_reads.fa.gz"]
    names += ["SIMULATED_0.5"]
    datas += ["/home/ashlemov/Git/ig_repertoire_constructor/test_dataset/igquast/SYNTHETIC_0.5/input_reads.fa.gz"]
    names += ["SYNTHETIC_0.5"]
    queue = zip(datas, names)

    def JOB(a):
        data, name = a
        hgc_estimator(data, name, tau=4)

    n_jobs = 1 if multiprocessing.cpu_count() <= 16 else 3
    Parallel(n_jobs=n_jobs)(delayed(JOB)(q) for q in queue)
    data = "/home/ashlemov/Git/ig_repertoire_constructor/test_dataset/igquast/age3/input_reads.fa.gz"
    name = "REAL_ALL"
    # hgc_estimator_all(data, name)
