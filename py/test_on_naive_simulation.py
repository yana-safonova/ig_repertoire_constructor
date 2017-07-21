#!/usr/bin/env python2

from simulate import *
from igquast_impl import get_clusters_sizes
import os.path
from joblib import Parallel, delayed
import multiprocessing
import tempfile

path_to_igquast = igrec_dir + "/igquast.py"


def run_and_quast_all(input_reads,
                      ideal_repertoire_fa,
                      ideal_repertoire_rcm,
                      out_dir,
                      threads=4,
                      rerun_mixcr=False,
                      rerun_igrec=False,
                      rerun_igquast=False,
                      do_not_run=False,
                      do_run_igrec_old=False,
                      do_run_igrec=True):
    import os.path
    import shutil

    mkdir_p(out_dir)

    class IgReCRun:

        def __init__(self,
                     name,
                     tau=4,
                     fillin=0.6,
                     min_sread_size=5,
                     max_votes=1000500,
                     trivial=False,
                     loci="all",
                     additional_args=" "):
            self.name = name
            self.tau = tau
            self.fillin = fillin
            self.trivial = trivial
            self.loci = loci
            self.min_sread_size = min_sread_size
            self.additional_args = additional_args
            self.max_votes = max_votes

        def run(self):
            additional_args = " --no-alignment --max-votes=%d" % self.max_votes
            if self.trivial:
                additional_args += " --create-triv-dec "
            additional_args += " " + self.additional_args
            run_igrec(input_reads,
                      tau=self.tau,
                      additional_args=additional_args,
                      min_fillin=self.fillin,
                      threads=threads,
                      loci=self.loci,
                      min_sread_size=self.min_sread_size,
                      output_dir=out_dir + "/" + self.name + "/")

    igrec_runs = []
    # igrec_runs.append(IgReCRun("igrec_trivial", trivial=True))
    igrec_runs.append(IgReCRun("igrec"))
    # igrec_runs.append(IgReCRun("igrec_tau3", tau=3))
    # igrec_runs.append(IgReCRun("igrec_nomsns", min_sread_size=1005000))
    # igrec_runs.append(IgReCRun("igrec_msns2", min_sread_size=2))
    # igrec_runs.append(IgReCRun("igrec_msns3", min_sread_size=3))
    # igrec_runs.append(IgReCRun("igrec_vote", max_votes=1))
    # igrec_runs.append(IgReCRun("igrec_vote2", max_votes=2))
    # igrec_runs.append(IgReCRun("igrec_trivial_tau3", tau=3, trivial=True))
    # igrec_runs.append(IgReCRun("igrec_trivial_tau2", tau=2, trivial=True))
    # igrec_runs.append(IgReCRun("igrec_trivial_tau1", tau=1, trivial=True))

    # igrec_runs.append(IgReCRun("igrec", additional_args="--debug"))
    # igrec_runs.append(IgReCRun("igrec_split", additional_args="--no-equal-compression --debug"))
    # igrec_runs.append(IgReCRun("igrec_split_tau3", tau=3, additional_args=" --no-equal-compression --debug"))
    # igrec_runs.append(IgReCRun("igrec_tau2", tau=2))
    # igrec_runs.append(IgReCRun("igrec_tau1", tau=1))

    # igrec_runs.append(IgReCRun("igrec_tau3_msns2", tau=3, min_sread_size=2))

    # igrec_runs.append(IgReCRun("igrec_tau3_vote", tau=3, max_votes=1))
    # igrec_runs.append(IgReCRun("igrec_tau3_vote2", tau=3, max_votes=2))
    # igrec_runs.append(IgReCRun("igrec_tau2_msns2", tau=2, min_sread_size=2))
    # igrec_runs.append(IgReCRun("igrec_tau1_msns2", tau=1, min_sread_size=2))

    # # igrec_runs.append(IgReCRun("igrec_msns3", min_sread_size=3))
    # # igrec_runs.append(IgReCRun("igrec_msns4", min_sread_size=4))
    # # igrec_runs.append(IgReCRun("igrec_f03", fillin=0.3))
    # # igrec_runs.append(IgReCRun("igrec_f075", fillin=0.75))
    # # igrec_runs.append(IgReCRun("igrec_f09", fillin=0.9))
    # # igrec_runs.append(IgReCRun("igrec_f095", fillin=0.95))
    # igrec_runs.append(IgReCRun("igrec_tau3", tau=3))
    # igrec_runs.append(IgReCRun("igrec_tau3_msns3", tau=3, min_sread_size=3))
    # igrec_runs.append(IgReCRun("igrec_tau3_f03", tau=3, fillin=0.3))
    # # igrec_runs.append(IgReCRun("igrec_tau3_f075", tau=3, fillin=0.75))
    # # igrec_runs.append(IgReCRun("igrec_tau3_f09", tau=3, fillin=0.9))
    # igrec_runs.append(IgReCRun("igrec_tau3_f095", tau=3, fillin=0.95))
    # igrec_runs.append(IgReCRun("igrec_tau2_msns2", tau=2, min_sread_size=2))
    # # igrec_runs.append(IgReCRun("igrec_tau2", tau=2))
    # # igrec_runs.append(IgReCRun("igrec_tau2_f03", tau=2, fillin=0.3))
    # # igrec_runs.append(IgReCRun("igrec_tau2_f075", tau=2, fillin=0.75))
    # # igrec_runs.append(IgReCRun("igrec_tau2_f09", tau=2, fillin=0.9))
    # igrec_runs.append(IgReCRun("igrec_tau2_f095", tau=2, fillin=0.95))
    # igrec_runs.append(IgReCRun("igrec_tau2_f099", tau=2, fillin=0.99))
    # # igrec_runs.append(IgReCRun("igrec_tau1", tau=1))
    # # igrec_runs.append(IgReCRun("igrec_tau1_f03", tau=1, fillin=0.3))
    # # igrec_runs.append(IgReCRun("igrec_tau1_f075", tau=1, fillin=0.75))
    # # igrec_runs.append(IgReCRun("igrec_tau1_f09", tau=1, fillin=0.9))
    # # igrec_runs.append(IgReCRun("igrec_tau1_f095", tau=1, fillin=0.95))

    if not do_not_run:
        if do_run_igrec:
            for run in igrec_runs:
                if rerun_igrec or not os.path.isfile(out_dir + "/" + run.name + "/final_repertoire.fa"):
                    run.run()

        if rerun_mixcr or not os.path.isfile(out_dir + "/mixcr2/final_repertoire.fa"):
            run_mixcr2(input_reads, threads=threads, output_dir=out_dir + "/mixcr2/", loci="all")

        if rerun_mixcr or not os.path.isfile(out_dir + "/mixcr2full/final_repertoire.fa"):
            run_mixcr2(input_reads, threads=threads, output_dir=out_dir + "/mixcr2full/", loci="all", region_from="FR1Begin", region_to="FR4End")

        if do_run_igrec_old:
            if rerun_igrec or not os.path.isfile(out_dir + "/" + "ig_repertoire_constructor" + "/final_repertoire.fa"):
                run_igrec_old(input_reads, threads=threads, output_dir=out_dir + "/ig_repertoire_constructor")
        # if rerun_mixcr or not os.path.isfile(out_dir + "/mixcr/final_repertoire.fa"):
        #     run_mixcr(input_reads, threads=threads, output_dir=out_dir + "/mixcr/", loci="all")

        # run_presto(input_reads, output_dir=out_dir + "/presto/")

        mkdir_p(out_dir + "/supernode")
        shutil.copy(out_dir + "/" + igrec_runs[0].name + "/supernode_repertoire.fa",
                    out_dir + "/supernode/final_repertoire.fa")
        shutil.copy(out_dir + "/" + igrec_runs[0].name + "/supernode_repertoire.rcm",
                    out_dir + "/supernode/final_repertoire.rcm")

    kinds = [run.name for run in igrec_runs] + ["supernode", "mixcr2", "mixcr2full"]

    if do_run_igrec_old:
        kinds += ["ig_repertoire_constructor"]

    for kind in kinds:
        args = {"ideal_repertoire_fa": ideal_repertoire_fa,
                "ideal_repertoire_rcm": ideal_repertoire_rcm,
                "input_reads": input_reads,
                "out_dir": out_dir,
                "kind": kind}
        if not rerun_igquast and os.path.isfile("%(out_dir)s/%(kind)s/aimquast/aimquast.json" % args):
            continue
        cmd = "gzip -f -k %(out_dir)s/%(kind)s/final_repertoire.fa" % args
        os.system(cmd)
        cmd = path_to_igquast + " -s %(input_reads)s \
            -r %(ideal_repertoire_fa)s -c %(out_dir)s/%(kind)s/final_repertoire.fa \
            -o %(out_dir)s/%(kind)s/aimquast --no-reference-free -F png,pdf,svg \
            --json %(out_dir)s/%(kind)s/aimquast/aimquast.json" % args

        # rcm = "%(out_dir)s/%(kind)s/final_repertoire.rcm" % args
        # if os.path.isfile(rcm):
        #     cmd += " -C %s" % rcm

        os.system(cmd)


if __name__ == "__main__":
    datasets = [igrec_dir + "/SIMULATED",
                igrec_dir + "/SYNTHETIC",
                igrec_dir + "/SIMTCR"]

    for dataset in datasets:
        mkdir_p(dataset)
        mkdir_p(dataset + "/data")

    if not os.path.isfile(datasets[0] + "/data/repertoire.fa.gz"):
        ig_simulator_output_dir = tempfile.mkdtemp()
        run_ig_simulator(ig_simulator_output_dir,
                         chain="HC", num_bases=1000, num_mutated=10000, repertoire_size=50000)
        fastx2fastx(ig_simulator_output_dir + "/ideal_repertoire.clusters.fa",
                    igrec_dir + "/SIMULATED/data/repertoire.fa.gz")
        rmdir(ig_simulator_output_dir)

    if not os.path.isfile(datasets[1] + "/data/repertoire.fa.gz"):
        try:
            convert_abvitro_to_repertoire("/Jake/data/input/ImmunoSeq/AbVitro/flu_time_course/FV/assembled_umis/21_assemble_combined.fastq",
                                          igrec_dir + "/SYNTHETIC/data/repertoire.fa.gz")
        except BaseException as ex:
            print ex
            print "Cannot convert FLU SYNTHETIC reperoire, file not found"

    if not os.path.isfile(datasets[2] + "/data/repertoire.fa.gz"):
        ig_simulator_output_dir = tempfile.mkdtemp()
        run_ig_simulator(ig_simulator_output_dir,
                         chain="HC", num_bases=150000, num_mutated=250000, repertoire_size=500000,
                         tcr=True)
        fastx2fastx(ig_simulator_output_dir + "/ideal_repertoire.clusters.fa",
                    igrec_dir + "/SIMTCR/data/repertoire.fa.gz")
        rmdir(ig_simulator_output_dir)


    for dataset in datasets:
        if not os.path.isfile(dataset + "/data/repertoire.fa.gz"):
            continue
        if not os.path.isfile(dataset + "/data/error_free_reads.fa.gz"):
            multiplex_repertoire(dataset + "/data/repertoire.fa.gz",
                                 dataset + "/data/error_free_reads.fa.gz")
            simulate_data_wo_errors(dataset + "/data/error_free_reads.fa.gz",
                                    dataset + "/data")

    print "Test datasets created!"

    lambdas = [0, 0.0625, 0.125, 0.25, 0.375, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
    for output_dir in datasets:
        print "Dataset:", output_dir
        if not os.path.isfile(output_dir + "/data/error_free_reads.fa.gz"):
            print "Not such dataset, skip!"
            continue
        is_synthetic = "SYNTHETIC" in output_dir
        is_simulated = not is_synthetic
        #
        min_error_interval = [0, 1]
        if "SYNTHETIC" not in output_dir:
            continue
        # if "SIMULATED" in output_dir or "SIMTCR" in output_dir:
        if "SIMULATED" in output_dir or "SYNTHETIC" in output_dir:
            seeds = range(10)
        else:
            seeds = [0]

        for min_error in min_error_interval:
            for error_rate in lambdas:
                for seed in seeds:
                    print "Seed", seed
                    out_dir = output_dir + "/errate_%0.4f" % error_rate if not min_error else output_dir + "/errate_%0.4f_woans" % error_rate
                    if len(seeds) != 1:
                        out_dir += "_seed_%d" % seed
                    if not os.path.isfile(out_dir + "/input_reads.fa.gz"):
                        # This code must not be run in parallel mode
                        simulate_data_from_dir(output_dir + "/data",
                                               out_dir,
                                               error_rate=error_rate,
                                               seed=seed,
                                               min_error=min_error,
                                               erroneous_site_len=300)

            def JOB(error_rate):
                print "Job:", error_rate
                for seed in seeds:
                    print "Seed", seed
                    out_dir = output_dir + "/errate_%0.4f" % error_rate if not min_error else output_dir + "/errate_%0.4f_woans" % error_rate
                    if len(seeds) != 1:
                        out_dir += "_seed_%d" % seed
                    if not os.path.isfile(out_dir + "/input_reads.fa.gz"):
                        simulate_data_from_dir(output_dir + "/data",
                                               out_dir,
                                               error_rate=error_rate,
                                               seed=seed,
                                               min_error=min_error,
                                               erroneous_site_len=300)

                    sizes = get_clusters_sizes(output_dir + "/data/repertoire.fa.gz")
                    print "Reference consists of %d clusters" % len(sizes)
                    print "Reference consists of %d large (>=5) clusters" % len([size for size in sizes if size >= 5])

                    run_and_quast_all(out_dir + "/input_reads.fa.gz",
                                      output_dir + "/data/repertoire.fa.gz",
                                      output_dir + "/data/repertoire.rcm", out_dir,
                                      rerun_mixcr=False,
                                      # do_run_igrec_old=True,
                                      do_run_igrec_old=False,
                                      threads=16)

            n_jobs = 1 if multiprocessing.cpu_count() <= 16 else 5

            _lambdas = lambdas if "TCR" not in output_dir else [0.5, 1, 2]
            Parallel(n_jobs=n_jobs)(delayed(JOB)(error_rate) for error_rate in _lambdas)
