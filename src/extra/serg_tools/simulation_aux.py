import os
import sys
from logging import INFO, ERROR
from string import join

current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = os.path.join(current_dir, os.pardir, os.pardir, os.pardir)
igrec_bin = os.path.join(igrec_dir, "build/release/bin")
migec_path = "/Marx/serg/soft/migec-1.2.4a"
mixcr_path = "/Marx/serg/soft/mixcr-2.0"

home_directory = os.path.abspath(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))))
pipeline_dir = os.path.join(home_directory, "py/pipeline/")
sys.path.append(pipeline_dir)
import support

sys.path.append(os.path.join(igrec_dir, "py"))
sys.path.append(os.path.join(igrec_dir, "py"))
from ash_python_utils import fastx2fastx
from simulate import run_mixcr2
from utils import fix_migec_mixcr_cluster_sizes


class RunParams:
    def __init__(self, data_path, pcr_error_rate, supernode_threshold, barcode_length):
        self.data_path = data_path
        self.pcr_error_rate = pcr_error_rate
        self.supernode_threshold = supernode_threshold
        self.barcode_length = barcode_length


class ShStep:
    def __init__(self, cwd, args):
        self.cwd = cwd
        self.cmdl = join(args, ' ')

    def Run(self, log):
        log.info("Running %s" % self.cmdl)
        try:
            support.sys_call(self.cmdl, log, self.cwd)
        except:
            log.error("Failed to run '%s':\n%s" % (self.cmdl, sys.exc_info()))
            return -1
        log.info("Returned 0")
        return 0


class PyStep:
    def __init__(self, desc, func):
        self.description = desc
        self.function = func

    def Run(self, log):
        log.info("Running python step: %s" % self.description)
        try:
            return_value = self.function()
            (log.info if return_value == 0 else log.error)("Returned %s" % return_value)
        except:
            log.error("Exception thrown")
            log.exception(sys.exc_info())
            return -1
        return return_value if return_value is not None else 0


def GetSimulateSteps(params, run_params):
    simulate_steps = [
        ShStep(None, ["%s/simulate_barcoded" % igrec_bin,
                      "--input-file %s/final_repertoire.fasta" % run_params.data_path,
                      "--output-dir %s/amplified" % run_params.data_path,
                      "--umi-length %d" % run_params.barcode_length,
                      "--pcr-error1 %f" % run_params.pcr_error_rate,
                      "--pcr-error2 %f" % run_params.pcr_error_rate,
                      # "--pcr-rate 0.001",
                      ]),
        ShStep(igrec_dir, ["%s/vj_finder" % igrec_bin,
                           "--input-file %s/amplified/repertoire_comp.fasta" % run_params.data_path,
                           "--output-dir %s/vjf_reference" % run_params.data_path,
                           "--loci IG",
                           "--threads %d" % params.threads
                           ])
    ]
    return simulate_steps


def GetBarigrecSteps(params, run_params):
    barigrec_steps = [
        ShStep(None, ["python -u %s/barcoded_igrec.py" % igrec_dir,
                      "-s %s/amplified/amplified.fasta" % run_params.data_path,
                      "--output %s/igrec_umi" % run_params.data_path,
                      "--loci IG",
                      "--threads %d" % params.threads,
                      "--igrec-tau 2",
                      "--min-super-read-size %d" % run_params.supernode_threshold,
                      "--no-compilation",
                      "--detect-chimeras",
                      "--clustering-thr 20"
                      ]),
        ShStep(None, ["python -u %s/py/drop_ns.py" % igrec_dir,
                      "-i %s/igrec_umi/final_repertoire/final_repertoire.fa" % run_params.data_path,
                      "-o %s/igrec_umi/final_repertoire.fa" % run_params.data_path,
                      "-r %s/igrec_umi/final_repertoire/final_repertoire.rcm" % run_params.data_path,
                      "-R %s/igrec_umi/final_repertoire.rcm" % run_params.data_path,
                      ]),
        ShStep(None, ["python -u %s/igquast.py" % igrec_dir,
                      "-s %s/amplified/amplified.fasta" % run_params.data_path,
                      "-c %s/igrec_umi/final_repertoire.fa" % run_params.data_path,
                      "-C %s/igrec_umi/final_repertoire.rcm" % run_params.data_path,
                      "-r %s/vjf_reference/cleaned_reads.fa" % run_params.data_path,
                      "-o %s/quast_barigrec" % run_params.data_path,
                      "--json %s/quast_barigrec/aimquast.json" % run_params.data_path,
                      "--reference-free",
                      "--rcm-based"
                      ])

        # ShStep(None, ["python -u %s/igquast.py" % igrec_dir,
        #  "-s %s/amplified/amplified.fasta" % simulation_params.data_path,
        #  "-c %s/igrec_umi/final_repertoire/final_repertoire.fa" % simulation_params.data_path,
        #  "-C %s/igrec_umi/final_repertoire/final_repertoire.rcm" % simulation_params.data_path,
        #  "-r %s/vjf_reference/cleaned_reads.fa" % simulation_params.data_path,
        #  "-o %s/quast_barigrec" % simulation_params.data_path,
        #  "--json %s/quast_barigrec/aimquast.json" % simulation_params.data_path,
        #  "--reference-free",
        #  "--rcm-based"
        #  ]),
    ]
    return barigrec_steps


def GetIgrecSteps(params, run_params):
    igrec_steps = [
        ShStep(None, ["python -u %s/igrec.py" % igrec_dir,
                      "-s %s/amplified/amplified.fasta" % run_params.data_path,
                      "-o %s/igrec" % run_params.data_path,
                      "--threads %d" % params.threads,
                      "--loci IGH",
                      "--debug"
                      ]),
        ShStep(None, ["python -u %s/py/drop_ns.py" % igrec_dir,
                      "-i %s/igrec/final_repertoire.fa" % run_params.data_path,
                      "-o %s/igrec/final_repertoire_non.fa" % run_params.data_path,
                      "-r %s/igrec/final_repertoire.rcm" % run_params.data_path,
                      "-R %s/igrec/final_repertoire_non.rcm" % run_params.data_path,
                      ]),
        ShStep(None, ["python -u %s/igquast.py" % igrec_dir,
                      "-s %s/amplified/amplified.fasta" % run_params.data_path,
                      "-c %s/igrec/final_repertoire_non.fa" % run_params.data_path,
                      "-C %s/igrec/final_repertoire_non.rcm" % run_params.data_path,
                      "-r %s/vjf_reference/cleaned_reads.fa" % run_params.data_path,
                      "-o %s/quast_igrec" % run_params.data_path,
                      "--json %s/quast_igrec/aimquast.json" % run_params.data_path,
                      "--reference-free",
                      "--rcm-based"
                      ]),

        # ShStep(None, ["python -u %s/igquast.py" %igrec_dir,
        #  "-s %s/amplified/amplified.fasta" % simulation_params.data_path,
        #  "-c %s/igrec/final_repertoire.fa" % simulation_params.data_path,
        #  "-C %s/igrec/final_repertoire.rcm" % simulation_params.data_path,
        #  "-r %s/vjf_reference/cleaned_reads.fa" % simulation_params.data_path,
        #  "-o %s/quast_igrec" % simulation_params.data_path,
        #  "--json %s/quast_igrec/aimquast.json" % simulation_params.data_path,
        #  "--reference-free",
        #  "--rcm-based"
        #  ]),
    ]
    return igrec_steps


def GetPrestoSteps(params, run_params):
    presto_steps = [
        ShStep(None, ["mkdir -p %s/presto" % run_params.data_path]),
        ShStep(None, ["python -u /Marx/ashlemov/Git/ig_repertoire_constructor/py/convertAGE2PRESTO.py",
                      "%s/amplified/amplified.fasta" % run_params.data_path,
                      "%s/presto/amplified_for_presto.fasta" % run_params.data_path
                      ]),
        ShStep("%s/presto" % run_params.data_path,
               [os.path.join(igrec_dir, "py/serg_tools/run_simple.sh"),
                "amplified_for_presto.fasta"
                ]),
        ShStep(None, ["python -u %s/py/convert_presto_to_quast.py" % igrec_dir,
                      "-r %s/presto/MS12_collapse-unique.fasta" % run_params.data_path,
                      "-o %s/presto/presto.fasta" % run_params.data_path
                      ]),
        ShStep(None, ["python -u %s/py/drop_ns.py" % igrec_dir,
                      "-i %s/presto/presto.fasta" % run_params.data_path,
                      "-o %s/presto/presto_non.fasta" % run_params.data_path
                      ]),
        ShStep(None, ["python -u %s/igquast.py" % igrec_dir,
                      "-s %s/amplified/amplified.fasta" % run_params.data_path,
                      "-c %s/presto/presto_non.fasta" % run_params.data_path,
                      "-r %s/vjf_reference/cleaned_reads.fa" % run_params.data_path,
                      "-o %s/quast_presto" % run_params.data_path,
                      "--json %s/quast_presto/aimquast.json" % run_params.data_path,
                      "--reference-free",
                      "--rcm-based"
                      ]),

        # ShStep(None, ["python -u %s/igquast.py" %igrec_dir,
        #  "-s %s/amplified/amplified.fasta" % simulation_params.data_path,
        #  "-c %s/presto/presto.fasta" % simulation_params.data_path,
        #  "-r %s/vjf_reference/cleaned_reads.fa" % simulation_params.data_path,
        #  "-o %s/quast_presto" % simulation_params.data_path,
        #  "--json %s/quast_presto/aimquast.json" % simulation_params.data_path,
        #  "--reference-free",
        #  "--rcm-based"
        #  ]),
    ]
    return presto_steps


def GetMigecSteps(params, run_params):
    migec_steps = [
        ShStep(igrec_dir, ["%s/vj_finder" % igrec_bin,
                           "--input-file %s/amplified/amplified.fasta" % run_params.data_path,
                           "--output-dir %s/vjf_amplified" % run_params.data_path,
                           "--loci IG",
                           "--threads %d" % params.threads
                           ]),
        PyStep("converting fasta file (%s) to fastq format (%s)" % (
            "%s/vjf_amplified/cleaned_reads.fa" % run_params.data_path,
            "%s/vjf_amplified/cleaned_reads.fastq" % run_params.data_path),
               lambda: fastx2fastx(
                   "%s/vjf_amplified/cleaned_reads.fa" % run_params.data_path,
                   "%s/vjf_amplified/cleaned_reads.fastq" % run_params.data_path,
                   50,
                   True)
               ),
        ShStep(None, ["python -u %s/py/convert_sim_to_migec.py" % igrec_dir,
                      "-r %s/vjf_amplified/cleaned_reads.fastq" % run_params.data_path,
                      "-o %s/vjf_amplified/migec.fastq" % run_params.data_path
                      ]),
        ShStep(None, ["java -jar %s/migec.jar Assemble" % migec_path,
                      "-c %s/vjf_amplified/migec.fastq" % run_params.data_path,
                      ".",
                      "%s/migec" % run_params.data_path
                      ]),
        PyStep("running MIXCR",
               lambda: run_mixcr2(
                   input_file="%s/migec/migec.t5.fastq.gz" % run_params.data_path,
                   output_dir="%s/migec/mixcr" % run_params.data_path,
                   threads=params.threads,
                   remove_tmp=False
               )),
        PyStep("fixing cluster sizes",
               lambda: fix_migec_mixcr_cluster_sizes(
                   input_file="%s/migec/mixcr/final_repertoire.fa" % run_params.data_path,
                   rcm_file="%s/migec/mixcr/final_repertoire.rcm" % run_params.data_path,
                   output_file="%s/migec/final_repertoire.fa" % run_params.data_path,
               )),

        # ShStep(None, ["java -jar %s/mixcr.jar" % mixcr_path,
        #         "align -p kaligner2",
        #         "--chains IGH",
        #         "%s/migec/migec.t5.fastq.gz" % simulation_params.data_path,
        #         "%s/migec/alignments.vdcja" % simulation_params.data_path
        #         ]),
        # ShStep(None, ["java -jar %s/mixcr.jar assemble" % mixcr_path,
        #         "-t %d" % threads,
        #         "-OassemblingFeatures=VDJRegion",
        #         "%s/migec/alignments.vdcja" % simulation_params.data_path,
        #         "%s/migec/clones.clns" % simulation_params.data_path
        #         ]),
        # ShStep(None, ["java -jar %s/mixcr.jar exportClones" % mixcr_path,
        #         "-t %d" % threads,
        #         "%s/migec/clones.clns" % simulation_params.data_path,
        #         "%s/migec/clones.txt" % simulation_params.data_path
        #         ]),
        # ShStep(None, ["python -u %s/py/convert_mixcr_to_quast.py" % igrec_dir,
        #         "-r %s/migec/clones.txt" % simulation_params.data_path,
        #         "-o %s/migec/clones.fasta" % simulation_params.data_path
        #         ]),

        # ShStep(None, ["gunzip",
        #  "--keep",
        #  "--force",
        #  "%s/migec/migec.t5.fastq.gz" % simulation_params.data_path,
        #  "> %s/migec/migec.t5.fastq" % simulation_params.data_path
        #  ]),
        # ShStep(None, ["sed",
        #  "'s/ /_/g'",
        #  "%s/migec/migec.t5.fastq" % simulation_params.data_path,
        #  "> %s/migec/migec.fastq" % simulation_params.data_path
        #  ]),
        # ShStep(None, ["python -u %s/py/convert_migec_to_trie.py" % igrec_dir,
        #  "-r %s/migec/migec.fastq" % simulation_params.data_path,
        #  "-o %s/migec/migec.fasta" % simulation_params.data_path
        #  ]),
        # ShStep(None, ["%s/py/ig_compress_equal_clusters.py" % igrec_dir,
        #  "%s/migec/migec.fasta" % simulation_params.data_path,
        #  "%s/migec/migec_compressed.fasta" % simulation_params.data_path,
        #  "--barcode"
        #  ]),

        # ShStep(None, ["python -u %s/igquast.py" % igrec_dir,
        #         "-s %s/amplified/amplified.fasta" % simulation_params.data_path,
        #         "-c %s/migec/clones.fasta" % simulation_params.data_path,
        #         "-r %s/vjf_reference/cleaned_reads.fa" % simulation_params.data_path,
        #         "-o %s/quast_migec" % simulation_params.data_path,
        #         "--json %s/quast_migec/aimquast.json" % simulation_params.data_path,
        #         "--reference-free",
        #         "--rcm-based"
        #         ]),
        ShStep(None, ["python -u %s/igquast.py" % igrec_dir,
                      "-s %s/amplified/amplified.fasta" % run_params.data_path,
                      "-c %s/migec/final_repertoire.fa" % run_params.data_path,
                      # "-C %s/migec/final_repertoire.rcm" % simulation_params.data_path,
                      "-r %s/vjf_reference/cleaned_reads.fa" % run_params.data_path,
                      "-o %s/quast_migec" % run_params.data_path,
                      "--json %s/quast_migec/aimquast.json" % run_params.data_path,
                      "--reference-free",
                      "--rcm-based"
                      ])
    ]
    return migec_steps
