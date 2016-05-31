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
from ash_python_utils import CreateLogger, AttachFileLogger, linear_search, idFormatByFileName, smart_open, md5_file, fq2fa



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



# def simulate_final_repertoire_fasta_to_merge_reads_final_and_rcm(input_file, output_file, error_rate):
if __name__ == "__main__":
    input_file = "/ssd/ig_simulator/ig_simulator_test/final_repertoire.fasta"
    simulated_repertoire_to_rcm(input_file, "/ssd/ig_repertoire_constructor/tmp_dir/ideal_final_repertoire.rcm")
    simulated_repertoire_to_final_repertoire(input_file, "/ssd/ig_repertoire_constructor/tmp_dir/ideal_final_repertoire.fa")
    jitt_fa_file(input_file, "/ssd/ig_repertoire_constructor/tmp_dir/merged_reads.fa", 2)
    jitt_fa_file(input_file, "/ssd/ig_repertoire_constructor/tmp_dir/merged_reads.fq", 2)
