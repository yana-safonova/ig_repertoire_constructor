import os
import shutil
import sys

from string import join

current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = os.path.join(current_dir, os.pardir, os.pardir, os.pardir)
sys.path.append(igrec_dir + "/src/extra/ash_python_utils/")

from ash_python_utils import fastx2fastx


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
        return_value = self.function()
        print "Returned", return_value
        return return_value if return_value is not None else 0


def run_sim_pipeline(data_path, pcr_error_rate, supernode_threshold, barcode_length, threads):
    home = "~/work" if os.path.exists("~/work") else "~"
    igrec = home + "/igrec"
    igrec_bin = igrec + "/build/release/bin"
    migec_path = "/Marx/serg/soft/migec-1.2.4a"
    mixcr_path = "/Marx/serg/soft/mixcr-2.0"

    simulate_steps = [
        ShStep(["%s/build/release/bin/simulate_barcoded" % igrec,
                "--input-file %s/final_repertoire.fasta" % data_path,
                "--output-file %s/amplified.fasta" % data_path,
                "--umi-length %d" % barcode_length,
                "--pcr-error1 %f" % pcr_error_rate,
                "--pcr-error2 %f" % pcr_error_rate,
                "--compressed-path %s/final_repertoire_comp.fasta" % data_path
                ]),
        ShStep(["cd %s &&" % igrec,
                "%s/vj_finder" % igrec_bin,
                "--input-file %s/final_repertoire_comp.fasta" % data_path,
                "--output-dir %s/vjf" % data_path,
                "--loci IG",
                "--threads %d" % threads
                ])
    ]

    barigrec_steps = [
        ShStep(["python %s/igrec_umi.py" % igrec,
                "-s %s/amplified.fasta" % data_path,
                "--output %s/igrec_umi" % data_path,
                "--loci IG",
                "--threads %d" % threads,
                "--igrec_tau 2",
                "--min-super-read-size %d" % supernode_threshold,
                "--no-compilation",
                "--detect-chimeras",
                "--clustering-thr 20"
                ]),
        ShStep(["python %s/py/drop_ns.py" % igrec,
                "-i %s/igrec_umi/final_repertoire/final_repertoire.fa" % data_path,
                "-o %s/igrec_umi/final_repertoire.fa" % data_path,
                "-r %s/igrec_umi/final_repertoire/final_repertoire.rcm" % data_path,
                "-R %s/igrec_umi/final_repertoire.rcm" % data_path,
                ]),
        ShStep(["python %s/aimquast.py" % igrec,
                "-s %s/amplified.fasta" % data_path,
                "-c %s/igrec_umi/final_repertoire.fa" % data_path,
                "-C %s/igrec_umi/final_repertoire.rcm" % data_path,
                "-r %s/vjf/cleaned_reads.fa" % data_path,
                "-o %s/quast" % data_path,
                "--reference-free",
                "--rcm-based"
                ])

        # ShStep(["python %s/aimquast.py" % igrec,
        #  "-s %s/amplified.fasta" % data_path,
        #  "-c %s/igrec_umi/final_repertoire/final_repertoire.fa" % data_path,
        #  "-C %s/igrec_umi/final_repertoire/final_repertoire.rcm" % data_path,
        #  "-r %s/vjf/cleaned_reads.fa" % data_path,
        #  "-o %s/quast" % data_path,
        #  "--reference-free",
        #  "--rcm-based"
        #  ]),
    ]

    igrec_steps = [
        ShStep(["python %s/igrec.py" % igrec,
                "-s %s/amplified.fasta" % data_path,
                "-o %s/igrec" % data_path,
                "--threads %d" % threads,
                "--loci IGH",
                "--debug"
                ]),
        ShStep(["python %s/py/drop_ns.py" % igrec,
                "-i %s/igrec/final_repertoire.fa" % data_path,
                "-o %s/igrec/final_repertoire_non.fa" % data_path,
                "-r %s/igrec/final_repertoire.rcm" % data_path,
                "-R %s/igrec/final_repertoire_non.rcm" % data_path,
                ]),
        ShStep(["python %s/aimquast.py" % igrec,
                "-s %s/amplified.fasta" % data_path,
                "-c %s/igrec/final_repertoire_non.fa" % data_path,
                "-C %s/igrec/final_repertoire_non.rcm" % data_path,
                "-r %s/vjf/cleaned_reads.fa" % data_path,
                "-o %s/quast_igrec" % data_path,
                "--reference-free",
                "--rcm-based"
                ]),

        # ShStep(["python %s/aimquast.py" %igrec,
        #  "-s %s/amplified.fasta" % data_path,
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
                "python /Marx/ashlemov/Git/ig_repertoire_constructor/py/convertAGE2PRESTO.py",
                "%s/amplified.fasta" % data_path,
                "%s/presto/amplified_for_presto.fasta" % data_path
                ]),
        ShStep(["cd %s/presto &&" % data_path,
                "../../run_simple.sh",
                "amplified_for_presto.fasta"
                ]),
        ShStep(["python %s/py/convert_presto_to_quast.py" % igrec,
                "-r %s/presto/MS12_collapse-unique.fasta" % data_path,
                "-o %s/presto/presto.fasta" % data_path
                ]),
        ShStep(["python %s/py/drop_ns.py" % igrec,
                "-i %s/presto/presto.fasta" % data_path,
                "-o %s/presto/presto_non.fasta" % data_path
                ]),
        ShStep(["python %s/aimquast.py" % igrec,
                "-s %s/amplified.fasta" % data_path,
                "-c %s/presto/presto_non.fasta" % data_path,
                "-r %s/vjf/cleaned_reads.fa" % data_path,
                "-o %s/quast_presto" % data_path,
                "--reference-free",
                "--rcm-based"
                ]),

        # ShStep(["python %s/aimquast.py" %igrec,
        #  "-s %s/amplified.fasta" % data_path,
        #  "-c %s/presto/presto.fasta" % data_path,
        #  "-r %s/vjf/cleaned_reads.fa" % data_path,
        #  "-o %s/quast_presto" % data_path,
        #  "--reference-free",
        #  "--rcm-based"
        #  ]),
    ]

    migec_steps = [
        ShStep(["cd %s &&" % igrec,
                "%s/vj_finder" % igrec_bin,
                "--input-file %s/amplified.fasta" % data_path,
                "--output-dir %s/vjf_input" % data_path,
                "--loci IG",
                "--threads %d" % threads
                ]),
        PyStep("converting fasta file (%s) to fastq format (%s)" % ("%s/vjf_input/cleaned_reads.fa" % data_path, "%s/vjf_input/cleaned_reads.fastq" % data_path),
               lambda:
               fastx2fastx(
                   "%s/vjf_input/cleaned_reads.fa" % data_path,
                   "%s/vjf_input/cleaned_reads.fastq" % data_path,
                   50,
                   True)
               ),
        ShStep(["python %s/py/convert_sim_to_migec.py" % igrec,
                "-r %s/vjf_input/cleaned_reads.fastq" % data_path,
                "-o %s/vjf_input/migec.fastq" % data_path
                ]),
        ShStep(["java -jar %s/migec.jar Assemble" % migec_path,
                "-c %s/vjf_input/migec.fastq" % data_path,
                ".",
                "%s/migec" % data_path
                ]),
        ShStep(["java -jar %s/mixcr.jar" % mixcr_path,
                "align -p kaligner2",
                "--chains IGH",
                "%s/migec/migec.t5.fastq.gz" % data_path,
                "%s/migec/alignments.vdcja" % data_path
                ]),
        ShStep(["java -jar %s/mixcr.jar assemble" % mixcr_path,
                "-t %d" % threads,
                "-OassemblingFeatures=VDJRegion",
                "%s/migec/alignments.vdcja" % data_path,
                "%s/migec/clones.clns" % data_path
                ]),
        ShStep(["java -jar %s/mixcr.jar exportClones" % mixcr_path,
                "-t %d" % threads,
                "%s/migec/clones.clns" % data_path,
                "%s/migec/clones.txt" % data_path
                ]),
        ShStep(["python %s/py/convert_mixcr_to_quast.py" % igrec,
                "-r %s/migec/clones.txt" % data_path,
                "-o %s/migec/clones.fasta" % data_path
                ]),

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
        # ShStep(["python %s/py/convert_migec_to_trie.py" % igrec,
        #  "-r %s/migec/migec.fastq" % data_path,
        #  "-o %s/migec/migec.fasta" % data_path
        #  ]),
        # ShStep(["%s/py/ig_compress_equal_clusters.py" % igrec,
        #  "%s/migec/migec.fasta" % data_path,
        #  "%s/migec/migec_compressed.fasta" % data_path,
        #  "--barcode"
        #  ]),
        ShStep(["python %s/aimquast.py" % igrec,
                "-s %s/amplified.fasta" % data_path,
                "-c %s/migec/clones.fasta" % data_path,
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
        if exit_status != 0:
            exit(exit_status)


def main():
    skip = 0 if len(sys.argv) < 3 else int(sys.argv[2])
    for supernode_threshold in [100000, 10, 5]:
        for barcode_length in [15, 9, 12]:
            for pcr_error_rate in [0.002, 0.0006, 0.0002]:
                if skip > 0:
                    skip -= 1
                    continue
                data_path = "%s/pcr_%f_super_%d_umi_%d" % (os.getcwd(), pcr_error_rate, supernode_threshold, barcode_length)
                if not os.path.exists(data_path):
                    os.makedirs(data_path)
                shutil.copyfile("final_repertoire.fasta", "%s/final_repertoire.fasta" % data_path)
                run_sim_pipeline(data_path, pcr_error_rate, supernode_threshold, barcode_length, int(sys.argv[1]))


if __name__ == '__main__':
    main()
