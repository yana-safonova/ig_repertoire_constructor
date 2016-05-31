#!/usr/bin/env python2

from argparse import ArgumentParser
from Bio import SeqIO
import sys
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import os
from matplotlib.backends.backend_pdf import PdfPages


current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = current_dir + "/../../../"
sys.path.append(igrec_dir + "/src/ig_tools/python_utils")
sys.path.append(igrec_dir + "/src/python_pipeline/")
from igblast_utils import ParseIgBlastOutput
import support
sys.path.append(igrec_dir + "/src/extra/ash_python_utils/")
from ash_python_utils import CreateLogger, AttachFileLogger, linear_search, idFormatByFileName, smart_open, md5_file, fq2fa, mkdir_p



def parse_final_repertoire_id(id):
    import re
    id = id.strip()

    m = re.match(r"^antibody_(\d+)_multiplicity_(\d+)_copy_(\d+)$", id)

    if m:
        g = m.groups()
        return int(g[0]), int(g[1]), int(g[2])
    else:
        return 1


assert parse_final_repertoire_id("antibody_1_multiplicity_1_copy_1") == (1, 1, 1)



def simulated_repertoire_to_rcm(input_file, rcm_file):
    with open(input_file) as fh, open(rcm_file, "w") as fout:
        for record in SeqIO.parse(fh, "fasta"):
            id = record.description
            cluster = str(parse_final_repertoire_id(id)[0])
            fout.write("%s\t%s\n" % (id, cluster))


def RC(l):
    import random

    S = set(list("ACTG"))
    s = S.difference([l])
    return random.choice(list(s))


def jitt_fa_file(input_file, output_file, error_rate=2, random_errors=True, min_error=0, erroneous_site_len=300, seed=None):
    import numpy as np
    from Bio import Seq
    import random

    output_format = idFormatByFileName(output_file)
    random.seed(seed)

    with smart_open(input_file) as fh, smart_open(output_file, "w") as fout:
        for record in SeqIO.parse(fh, "fasta"):
            n_errors = np.random.poisson(error_rate, 1)[0] if random_errors else error_rate
            if n_errors < min_error:
                n_errors = min_error

            positions = random.sample(range(min(len(record.seq), erroneous_site_len)), n_errors)
            s = list(str(record.seq))
            for pos in positions:
                s[pos] = RC(s[pos])
            record.letter_annotations = {}
            record.seq = Seq.Seq("".join(s))

            if output_format == "fastq":
                # record.letter_annotations["phred_quality"] = [random.randint(30, 50) for _ in xrange(len(record))]  # TODO Check it out
                record.letter_annotations["phred_quality"] = [50] * len(record)  # TODO Check it out

            SeqIO.write(record, fout, output_format)


def simulated_repertoire_to_final_repertoire(input_file, output_file):
    import random

    output_format = idFormatByFileName(output_file)

    with smart_open(input_file) as fh, smart_open(output_file, "w") as fout:
        for record in SeqIO.parse(fh, "fasta"):
            id = record.description
            cluster, size, copy = parse_final_repertoire_id(id)
            if copy == 1:
                record.id = record.description = "cluster___%d___size___%d" % (cluster, size)
                record.letter_annotations = {}

                if output_format == "fastq":
                    # record.letter_annotations["phred_quality"] = [random.randint(30, 50) for _ in xrange(len(record))]  # TODO Check it out
                    record.letter_annotations["phred_quality"] = [50] * len(record)

                SeqIO.write(record, fout, output_format)


def simulate_data(input_file, output_dir, **kwargs):
    simulated_repertoire_to_rcm(input_file, "%s/ideal_final_repertoire.rcm" % output_dir)
    simulated_repertoire_to_final_repertoire(input_file, "%s/ideal_final_repertoire.fa" % output_dir)
    jitt_fa_file(input_file, "%s/merged_reads.fq" % output_dir, **kwargs)


path_to_ig_simulator = "/ssd/ig_simulator/"
path_to_mixcr = "TBD"
path_to_igrec = igrec_dir


def run_ig_simulator(output_dir, log, chain="HC", num_bases=100, num_mutated=1000, reprtoire_size=5000):
    assert chain in ["HC", "LC"]

    args = {"path": path_to_ig_simulator,
            "output_dir": output_dir,
            "chain": chain,
            "num_bases": num_bases,
            "num_mutated": num_mutated,
            "reprtoire_size": reprtoire_size}
    support.sys_call("%{path}s/ig_simulator.py --chain-type %{chain}s --num-bases %{num_bases}d --num-mutated %{num_mutated}d --repertoire-size %{reprtoire_size}d -o %s{output_dir}s --skip-drawing" % args,
                     log=log)


def convert_mixcr_output_to_igrec(input_file, output_file):
    with smart_open(input_file) as fh, smart_open(output_file, "w") as fout:
        # Skip header
        fh.next()

        for i, line in enumerate(fh):
            seq, size = line.strip().split()
            size = int(size)
            fout.write(">cluster___%d___size___%d\n" % (i, size))
            fout.write(seq + "\n")


def run_igreg(input_file, output_dir, log, tau=4, loci="IGH"):
    args = {"path": path_to_igrec,
            "tau": tau,
            "loci": loci,
            "input_file": input_file,
            "output_dir": output_dir}
    support.sys_call("%{path}s/igrec.py --tau=%{tau}d -t 36 --loci %{loci}s -s %{input_file}s -o %{output_dir}s" % args)


def fa2fq(input_file, output_file):
    import random

    with smart_open(input_file) as fh, smart_open(output_file, "w") as fout:
        for record in SeqIO.parse(fh, "fasta"):
            record.letter_annotations["phred_quality"] = [50] * len(record)
            # record.letter_annotations["phred_quality"] = [random.randint(30, 50) for _ in xrange(len(record))]
            SeqIO.write(record, fout, "fastq")


def run_mixcr(input_file, output_dir, log, loci="IGH"):
    mkdir_p(output_dir)

    if idFormatByFileName(input_file) == "fasta":
        input_file_fq = "%s/input_reads.fq" % output_dir
        fa2fq(input_file, input_file_fq)
        input_file = input_file_fq


    args = {"path": path_to_mixcr,
            "mixcr_cmd": "java -jar %s/mixcr.jar" % path_to_mixcr,
            "output_dir": output_dir}

    support.sys_call("%{mixcr_cmd}s  align -f -g -r %{output_dir}s/align_report.txt --loci %{loci}s --noMerge -OvParameters.geneFeatureToAlign=VTranscript %{input_file}s %{output_dir}s/mixcr.vdjca" % args,
                     log=log)
    support.sys_call("%{mixcr_cmd}s assemble -f -r %{output_dir}s/assemble_report.txt -OassemblingFeatures=\"{FR1Begin:CDR3End}\" %{output_dir}s/mixcr.vdjca %{output_dir}s/mixcr.clns" % args,
                     log=log)
    os.system("%{mixcr_cmd}s exportClones -pf %{path}s/preset.pf -f --no-spaces %{output_dir}s/mixcr.clns %{output_dir}s/mixcr.txt" % args)

    convert_mixcr_output_to_igrec("%{output_dir}s/mixcr.txt" % args, "%{output_dir}s/mixcr_final.fa" % args)


if __name__ == "__main__":
    log = CreateLogger("aimQUAST")

    input_file = "/ssd/ig_simulator/ig_simulator_test/final_repertoire.fasta"
    simulate_data(input_file, "/ssd/ig_repertoire_constructor/tmp_dir", error_rate=2)
