import os
import sys

import logging

word_dir = os.getcwd()
igrec_dir = os.path.join(word_dir, 'env')
igrec_bin = os.path.join(igrec_dir, 'build/release/bin')

from simulation_aux import ShStep, PyStep, migec_path

sys.path.append(os.path.join(igrec_dir, "py"))
sys.path.append(os.path.join(igrec_dir, "py"))
from ash_python_utils import fastx2fastx
from simulate import run_mixcr2
from utils import fix_migec_mixcr_cluster_sizes


threads = 32

steps = [
    ShStep(igrec_dir, ["%s/vj_finder" % igrec_bin,
                       "--input-file %s" % "/Marx/serg/data/age/3/age_ig_s3_R12.fastq",
                       "--output-dir %s/vjf_amplified" % word_dir,
                       "--loci IG",
                       "--threads %d" % threads
                       ]),
    PyStep("converting fasta file (%s) to fastq format (%s)" % (
        "%s/vjf_amplified/cleaned_reads.fa" % word_dir,
        "%s/vjf_amplified/cleaned_reads.fastq" % word_dir),
           lambda: fastx2fastx(
               "%s/vjf_amplified/cleaned_reads.fa" % word_dir,
               "%s/vjf_amplified/cleaned_reads.fastq" % word_dir,
               50,
               True)
           ),
    ShStep(None, ["python -u %s/py/convert_sim_to_migec.py" % igrec_dir,
                  "-r %s/vjf_amplified/cleaned_reads.fastq" % word_dir,
                  "-o %s/vjf_amplified/migec.fastq" % word_dir
                  ]),
    ShStep(None, ["java -jar %s/migec.jar Assemble" % migec_path,
                  "-c %s/vjf_amplified/migec.fastq" % word_dir,
                  ".",
                  "%s/migec" % word_dir
                  ]),
    PyStep("running MIXCR",
           lambda: run_mixcr2(
               input_file="%s/migec/migec.t5.fastq.gz" % word_dir,
               output_dir="%s/migec/mixcr" % word_dir,
               threads=threads,
               remove_tmp=False
           )),
    PyStep("fixing cluster sizes",
           lambda: fix_migec_mixcr_cluster_sizes(
               input_file="%s/migec/mixcr/final_repertoire.fa" % word_dir,
               rcm_file="%s/migec/mixcr/final_repertoire.rcm" % word_dir,
               output_file="%s/migec/final_repertoire.fa" % word_dir,
           ))
]


def CreateLogger(output_dir):
    log = logging.getLogger('run_sim')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    log_filename = os.path.join(output_dir, "sim.log")
    log_handler = logging.FileHandler(log_filename, mode='a' if os.path.exists(log_filename) else 'w')
    log.addHandler(log_handler)
    log.info("Log will be written to " + log_filename + "\n")
    return log


log = CreateLogger(word_dir)

for step in steps:
    exit_status = step.Run(log)
    if exit_status != 0:
        exit(exit_status)
