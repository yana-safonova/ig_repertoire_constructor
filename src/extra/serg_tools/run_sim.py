#!/usr/bin/python -u

import os
import shutil
import sys

from string import join

from py.plot_umi import plot_sens_prec_umi

current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = os.path.join(current_dir, os.pardir, os.pardir, os.pardir)
sys.path.append(igrec_dir + "/src/extra/ash_python_utils/")
sys.path.append(igrec_dir + "/py/")

from ash_python_utils import fastx2fastx
from simulate import run_mixcr2
from utils import fix_migec_mixcr_cluster_sizes


class ShStep:
    def __init__(self, args):
        self.cmdl = join(args, ' ')

    def Run(self):
        print "Running %s" % self.cmdl
        exit_status = os.system(self.cmdl)
        print "Returned", exit_status
        return exit_status


class PyStep:
    def __init__(self, desc, func):
        self.description = desc
        self.function = func

    def Run(self):
        print "Running python step: %s" % self.description
        try:
            return_value = self.function()
            print "Returned", return_value
        except:
            return -1
        return return_value if return_value is not None else 0


def run_sim_pipeline(data_path, pcr_error_rate, supernode_threshold, barcode_length, threads, exit_on_error):
    igrec_bin = igrec_dir + "/build/release/bin"
    migec_path = "/Marx/serg/soft/migec-1.2.4a"
    mixcr_path = "/Marx/serg/soft/mixcr-2.0"

    simulate_steps = [
        ShStep(["%s/simulate_barcoded" % igrec_bin,
                "--input-file %s/final_repertoire.fasta" % data_path,
                "--output-dir %s/amplified" % data_path,
                "--umi-length %d" % barcode_length,
                "--pcr-error1 %f" % pcr_error_rate,
                "--pcr-error2 %f" % pcr_error_rate,
                # "--pcr-rate 0.001",
                ]),
        ShStep(["cd %s &&" % igrec_dir,
                "%s/vj_finder" % igrec_bin,
                "--input-file %s/amplified/repertoire_comp.fasta" % data_path,
                "--output-dir %s/vjf" % data_path,
                "--loci IG",
                "--threads %d" % threads
                ])
    ]

    barigrec_steps = [
        ShStep(["python -u %s/igrec_umi.py" % igrec_dir,
                "-s %s/amplified/amplified.fasta" % data_path,
                "--output %s/igrec_umi" % data_path,
                "--loci IG",
                "--threads %d" % threads,
                "--igrec_tau 2",
                "--min-super-read-size %d" % supernode_threshold,
                "--no-compilation",
                "--detect-chimeras",
                "--clustering-thr 20"
                ]),
        ShStep(["python -u %s/py/drop_ns.py" % igrec_dir,
                "-i %s/igrec_umi/final_repertoire/final_repertoire.fa" % data_path,
                "-o %s/igrec_umi/final_repertoire.fa" % data_path,
                "-r %s/igrec_umi/final_repertoire/final_repertoire.rcm" % data_path,
                "-R %s/igrec_umi/final_repertoire.rcm" % data_path,
                ]),
        ShStep(["python -u %s/aimquast.py" % igrec_dir,
                "-s %s/amplified/amplified.fasta" % data_path,
                "-c %s/igrec_umi/final_repertoire.fa" % data_path,
                "-C %s/igrec_umi/final_repertoire.rcm" % data_path,
                "-r %s/vjf/cleaned_reads.fa" % data_path,
                "-o %s/quast" % data_path,
                "--reference-free",
                "--rcm-based"
                ])

        # ShStep(["python -u %s/aimquast.py" % igrec_dir,
        #  "-s %s/amplified/amplified.fasta" % data_path,
        #  "-c %s/igrec_umi/final_repertoire/final_repertoire.fa" % data_path,
        #  "-C %s/igrec_umi/final_repertoire/final_repertoire.rcm" % data_path,
        #  "-r %s/vjf/cleaned_reads.fa" % data_path,
        #  "-o %s/quast" % data_path,
        #  "--reference-free",
        #  "--rcm-based"
        #  ]),
    ]

    igrec_steps = [
        ShStep(["python -u %s/igrec.py" % igrec_dir,
                "-s %s/amplified/amplified.fasta" % data_path,
                "-o %s/igrec" % data_path,
                "--threads %d" % threads,
                "--loci IGH",
                "--debug"
                ]),
        ShStep(["python -u %s/py/drop_ns.py" % igrec_dir,
                "-i %s/igrec/final_repertoire.fa" % data_path,
                "-o %s/igrec/final_repertoire_non.fa" % data_path,
                "-r %s/igrec/final_repertoire.rcm" % data_path,
                "-R %s/igrec/final_repertoire_non.rcm" % data_path,
                ]),
        ShStep(["python -u %s/aimquast.py" % igrec_dir,
                "-s %s/amplified/amplified.fasta" % data_path,
                "-c %s/igrec/final_repertoire_non.fa" % data_path,
                "-C %s/igrec/final_repertoire_non.rcm" % data_path,
                "-r %s/vjf/cleaned_reads.fa" % data_path,
                "-o %s/quast_igrec" % data_path,
                "--reference-free",
                "--rcm-based"
                ]),

        # ShStep(["python -u %s/aimquast.py" %igrec_dir,
        #  "-s %s/amplified/amplified.fasta" % data_path,
        #  "-c %s/igrec/final_repertoire.fa" % data_path,
        #  "-C %s/igrec/final_repertoire.rcm" % data_path,
        #  "-r %s/vjf/cleaned_reads.fa" % data_path,
        #  "-o %s/quast_igrec" % data_path,
        #  "--reference-free",
        #  "--rcm-based"
        #  ]),
    ]

    presto_steps = [
        ShStep(["mkdir -p %s/presto &&" % data_path,
                "python -u /Marx/ashlemov/Git/ig_repertoire_constructor/py/convertAGE2PRESTO.py",
                "%s/amplified/amplified.fasta" % data_path,
                "%s/presto/amplified_for_presto.fasta" % data_path
                ]),
        ShStep(["cd %s/presto &&" % data_path,
                "../../run_simple.sh",
                "amplified_for_presto.fasta"
                ]),
        ShStep(["python -u %s/py/convert_presto_to_quast.py" % igrec_dir,
                "-r %s/presto/MS12_collapse-unique.fasta" % data_path,
                "-o %s/presto/presto.fasta" % data_path
                ]),
        ShStep(["python -u %s/py/drop_ns.py" % igrec_dir,
                "-i %s/presto/presto.fasta" % data_path,
                "-o %s/presto/presto_non.fasta" % data_path
                ]),
        ShStep(["python -u %s/aimquast.py" % igrec_dir,
                "-s %s/amplified/amplified.fasta" % data_path,
                "-c %s/presto/presto_non.fasta" % data_path,
                "-r %s/vjf/cleaned_reads.fa" % data_path,
                "-o %s/quast_presto" % data_path,
                "--reference-free",
                "--rcm-based"
                ]),

        # ShStep(["python -u %s/aimquast.py" %igrec_dir,
        #  "-s %s/amplified/amplified.fasta" % data_path,
        #  "-c %s/presto/presto.fasta" % data_path,
        #  "-r %s/vjf/cleaned_reads.fa" % data_path,
        #  "-o %s/quast_presto" % data_path,
        #  "--reference-free",
        #  "--rcm-based"
        #  ]),
    ]

    migec_steps = [
        ShStep(["cd %s &&" % igrec_dir,
                "%s/vj_finder" % igrec_bin,
                "--input-file %s/amplified/amplified.fasta" % data_path,
                "--output-dir %s/vjf_input" % data_path,
                "--loci IG",
                "--threads %d" % threads
                ]),
        PyStep("converting fasta file (%s) to fastq format (%s)" % ("%s/vjf_input/cleaned_reads.fa" % data_path, "%s/vjf_input/cleaned_reads.fastq" % data_path),
               lambda: fastx2fastx(
                   "%s/vjf_input/cleaned_reads.fa" % data_path,
                   "%s/vjf_input/cleaned_reads.fastq" % data_path,
                   50,
                   True)
               ),
        ShStep(["python -u %s/py/convert_sim_to_migec.py" % igrec_dir,
                "-r %s/vjf_input/cleaned_reads.fastq" % data_path,
                "-o %s/vjf_input/migec.fastq" % data_path
                ]),
        ShStep(["java -jar %s/migec.jar Assemble" % migec_path,
                "-c %s/vjf_input/migec.fastq" % data_path,
                ".",
                "%s/migec" % data_path
                ]),
        PyStep("running MIXCR",
               lambda: run_mixcr2(
                   input_file = "%s/migec/migec.t5.fastq.gz" % data_path,
                   output_dir = "%s/migec/mixcr" % data_path,
                   threads = threads,
                   remove_tmp = False
               )),
        PyStep("fixing cluster sizes",
               lambda: fix_migec_mixcr_cluster_sizes(
                   input_file = "%s/migec/mixcr/final_repertoire.fa" % data_path,
                   rcm_file = "%s/migec/mixcr/final_repertoire.rcm" % data_path,
                   output_file = "%s/migec/final_repertoire.fa" % data_path,
               )),

        # ShStep(["java -jar %s/mixcr.jar" % mixcr_path,
        #         "align -p kaligner2",
        #         "--chains IGH",
        #         "%s/migec/migec.t5.fastq.gz" % data_path,
        #         "%s/migec/alignments.vdcja" % data_path
        #         ]),
        # ShStep(["java -jar %s/mixcr.jar assemble" % mixcr_path,
        #         "-t %d" % threads,
        #         "-OassemblingFeatures=VDJRegion",
        #         "%s/migec/alignments.vdcja" % data_path,
        #         "%s/migec/clones.clns" % data_path
        #         ]),
        # ShStep(["java -jar %s/mixcr.jar exportClones" % mixcr_path,
        #         "-t %d" % threads,
        #         "%s/migec/clones.clns" % data_path,
        #         "%s/migec/clones.txt" % data_path
        #         ]),
        # ShStep(["python -u %s/py/convert_mixcr_to_quast.py" % igrec_dir,
        #         "-r %s/migec/clones.txt" % data_path,
        #         "-o %s/migec/clones.fasta" % data_path
        #         ]),

        # ShStep(["gunzip",
        #  "--keep",
        #  "--force",
        #  "%s/migec/migec.t5.fastq.gz" % data_path,
        #  "> %s/migec/migec.t5.fastq" % data_path
        #  ]),
        # ShStep(["sed",
        #  "'s/ /_/g'",
        #  "%s/migec/migec.t5.fastq" % data_path,
        #  "> %s/migec/migec.fastq" % data_path
        #  ]),
        # ShStep(["python -u %s/py/convert_migec_to_trie.py" % igrec_dir,
        #  "-r %s/migec/migec.fastq" % data_path,
        #  "-o %s/migec/migec.fasta" % data_path
        #  ]),
        # ShStep(["%s/py/ig_compress_equal_clusters.py" % igrec_dir,
        #  "%s/migec/migec.fasta" % data_path,
        #  "%s/migec/migec_compressed.fasta" % data_path,
        #  "--barcode"
        #  ]),

        # ShStep(["python -u %s/aimquast.py" % igrec_dir,
        #         "-s %s/amplified/amplified.fasta" % data_path,
        #         "-c %s/migec/clones.fasta" % data_path,
        #         "-r %s/vjf/cleaned_reads.fa" % data_path,
        #         "-o %s/quast_migec" % data_path,
        #         "--reference-free",
        #         "--rcm-based"
        #         ]),
        ShStep(["python -u %s/aimquast.py" % igrec_dir,
                "-s %s/amplified/amplified.fasta" % data_path,
                "-c %s/migec/final_repertoire.fa" % data_path,
                # "-C %s/migec/final_repertoire.rcm" % data_path,
                "-r %s/vjf/cleaned_reads.fa" % data_path,
                "-o %s/quast_migec" % data_path,
                "--reference-free",
                "--rcm-based"
                ])
    ]

    steps = []
    steps.extend(simulate_steps)
    steps.extend(barigrec_steps)
    steps.extend(igrec_steps)
    steps.extend(presto_steps)
    steps.extend(migec_steps)

    for step in steps:
        exit_status = step.Run()
        if exit_status != 0 and exit_on_error:
            exit(exit_status)


def main():
    print "Starting run_sum ", sys.argv

    skip = 0 if len(sys.argv) < 3 else int(sys.argv[2])
    exit_on_error = 1 if len(sys.argv) < 4 else int(sys.argv[3])
    base_data_path = os.getcwd()
    # for supernode_threshold in [100000, 10, 5]:
    for supernode_threshold in [100000]:
        # for barcode_length in [15, 9, 12]:
        for barcode_length in [15]:
            for pcr_error_rate in [0.006, 0.0006, 0.0025]:
                if skip > 0:
                    skip -= 1
                    continue
                data_path = "%s/pcr_%f_super_%d_umi_%d" % (base_data_path, pcr_error_rate, supernode_threshold, barcode_length)
                if not os.path.exists(data_path):
                    os.makedirs(data_path)
                shutil.copyfile("final_repertoire.fasta", "%s/final_repertoire.fasta" % data_path)
                run_sim_pipeline(data_path, pcr_error_rate, supernode_threshold, barcode_length, int(sys.argv[1]), exit_on_error)

    PyStep("Drawing plots",
           lambda: plot_sens_prec_umi(base_data_path))


if __name__ == '__main__':
    main()
