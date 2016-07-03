#!/usr/bin/env python2

from simulate import *
from aimquast_impl import get_clusters_sizes
import os


path_to_aimquast = igrec_dir + "/aimquast.py"


def run_and_quast_all(input_reads,
                      ideal_repertoire_fa,
                      ideal_repertoire_rcm,
                      out_dir,
                      threads=4,
                      rerun_mixcr=False,
                      do_not_run=False):
    import os.path

    mkdir_p(out_dir)

    if not do_not_run:
        run_igrec(input_reads,
                  threads=threads,
                  output_dir=out_dir + "/igrec/")

        run_igrec(input_reads,
                  additional_args=" --create-triv-dec",
                  threads=threads,
                  output_dir=out_dir + "/igrec_trivial/")

        run_igrec(input_reads,
                  min_fillin=0.3,
                  threads=threads,
                  output_dir=out_dir + "/igrec_f03/")

        run_igrec(input_reads,
                  tau=3,
                  threads=threads,
                  output_dir=out_dir + "/igrec_tau3/")

        run_igrec(input_reads,
                  tau=3,
                  min_fillin=0.3,
                  threads=threads,
                  output_dir=out_dir + "/igrec_tau3_f03/")

        run_igrec(input_reads,
                  tau=3,
                  min_fillin=0.75,
                  threads=threads,
                  output_dir=out_dir + "/igrec_tau3_f075/")

        if rerun_mixcr or not os.path.isfile(out_dir + "/mixcr/final_repertoire.fa"):
            run_mixcr(input_reads, threads=threads, output_dir=out_dir + "/mixcr/")

        # run_igrec(input_reads,
        #           additional_args=" --create-triv-dec",
        #           threads=threads,
        #           output_dir=out_dir + "/igrec_trivial/")
        #
        # run_igrec(input_reads,
        #           min_fillin=0.3,
        #           threads=threads,
        #           output_dir=out_dir + "/igrec_f03/")
        #
        # run_igrec(input_reads,
        #           tau=3,
        #           threads=threads,
        #           output_dir=out_dir + "/igrec_tau3/")
        #
        # run_igrec(input_reads,
        #           tau=3,
        #           min_fillin=0.3,
        #           threads=threads,
        #           output_dir=out_dir + "/igrec_tau3_f03/")
        #
        # run_igrec(input_reads,
        #           tau=3,
        #           min_fillin=0.75",
        #           threads=threads,
        #           output_dir=out_dir + "/igrec_tau3_f075/")
        #
        # if rerun_mixcr or not os.path.isfile(out_dir + "/mixcr/final_repertoire.fa"):
        #     run_mixcr(input_reads, threads=threads, output_dir=out_dir + "/mixcr/")

        run_igrec(input_reads,
                  tau=3,
                  min_fillin=0.9,
                  threads=threads,
                  output_dir=out_dir + "/igrec_tau3_f09/")

        run_igrec(input_reads,
                  tau=3,
                  min_fillin=0.95,
                  threads=threads,
                  output_dir=out_dir + "/igrec_tau3_f095/")

        run_igrec(input_reads,
                  tau=2,
                  min_fillin=0.6,
                  threads=threads,
                  output_dir=out_dir + "/igrec_tau2_f06/")

        run_igrec(input_reads,
                  tau=2,
                  min_fillin=0.9,
                  threads=threads,
                  output_dir=out_dir + "/igrec_tau2_f09/")

        run_igrec(input_reads,
                  tau=2,
                  min_fillin=0.3,
                  threads=threads,
                  output_dir=out_dir + "/igrec_tau2_f03/")

        run_igrec(input_reads,
                  tau=1,
                  min_fillin=0.6,
                  threads=threads,
                  output_dir=out_dir + "/igrec_tau1_f06/")

    kinds = ["igrec", "igrec_trivial", "igrec_f03", "igrec_tau3", "igrec_tau3_f03", "igrec_tau3_f075", "mixcr"]
    kinds2 = ["igrec_tau3_f09", "igrec_tau3_f095", "igrec_tau2_f06", "igrec_tau2_f09", "igrec_tau2_f03", "igrec_tau1_f06"]

    for kind in kinds + kinds2:
        args = {"ideal_repertoire_fa": ideal_repertoire_fa,
                "ideal_repertoire_rcm": ideal_repertoire_rcm,
                "input_reads": input_reads,
                "out_dir": out_dir,
                "kind": kind}
        cmd = path_to_aimquast + " -s %(input_reads)s -r %(ideal_repertoire_fa)s -R %(ideal_repertoire_rcm)s -c %(out_dir)s/%(kind)s/final_repertoire.fa -o %(out_dir)s/%(kind)s/aimquast --no-reference-free -F png,pdf" % args

        rcm = "%(out_dir)s/%(kind)s/final_repertoire.rcm" % args
        if os.path.isfile(rcm):
            cmd += " -C %s" % rcm

        os.system(cmd)


def run_and_quast_all_three(input_dir, threads=16):
    # run_and_quast_all(input_dir + "/input1.fa.gz",
    #                   input_dir + "/repertoire1.fa.gz",
    #                   input_dir + "/repertoire1.rcm",
    #                   out_dir=input_dir + "/filtering1",
    #                   threads=threads)
    #
    # run_and_quast_all(input_dir + "/input2.fa.gz",
    #                   input_dir + "/repertoire2.fa.gz",
    #                   input_dir + "/repertoire2.rcm",
    #                   out_dir=input_dir + "/filtering2",
    #                   threads=threads)

    run_and_quast_all(input_dir + "/input3.fa.gz",
                      input_dir + "/repertoire3.fa.gz",
                      input_dir + "/repertoire3.rcm",
                      out_dir=input_dir + "/filtering3",
                      threads=threads)


if __name__ == "__main__":
    import os.path

    datasets = ["AGE1", "AGE2", "AGE3", "AGE4", "AGE5", "AGE6", "AGE7", "AGE8", "AGE9", "MG91M", "HD09M"]
    datasets = ["AGE3", "MG91M", "HD09M", "AGE7"]
    # datasets = ["AGE3", "MG91M", "HD09M"]
    datasets = ["AGE3"]
    datasets = []
    datasets = ["FLU_FV_27"]
    datasets = ["FLU_FV_21", "FLU_FV_22", "FLU_FV_23", "FLU_FV_27"]
    for dataset in datasets:
        ddir = igrec_dir + "/src/extra/ig_quast_tool/" + dataset
        if os.path.isfile(ddir + "/input1.fa.gz"):
            run_and_quast_all_three(ddir)

    ig_simulator_output_dir = "/tmp/ig_simulator"
    output_dir = igrec_dir + "/various_error_rate"
    run_ig_simulator(ig_simulator_output_dir,
                     chain="HC", num_bases=100, num_mutated=1000, reprtoire_size=5000)

    lambdas = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4]
    # lambdas = [0, 1, 2]
    for min_error in [0, 1]:
        for error_rate in lambdas:
            out_dir = output_dir + "/errate_%0.2f" % error_rate if not min_error else output_dir + "/errate_%0.2f_woans" % error_rate
            simulate_data(ig_simulator_output_dir + "/final_repertoire.fasta",
                          out_dir,
                          error_rate=error_rate,
                          seed=0,
                          min_error=min_error,
                          erroneous_site_len=300)

            sizes = get_clusters_sizes(out_dir + "/ideal_final_repertoire.fa")
            print "Reference consists of %d clusters" % len(sizes)
            print "Reference consists of %d large (>=5) clusters" % len([size for size in sizes if size >= 5])
            sys.exit()

            run_and_quast_all(out_dir + "/merged_reads.fa",
                              out_dir + "/ideal_final_repertoire.fa",
                              out_dir + "/ideal_final_repertoire.rcm", out_dir,
                              rerun_mixcr=True)
