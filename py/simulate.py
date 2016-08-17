#!/usr/bin/env python2

import os
import sys
from Bio import SeqIO

current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = current_dir + "/../"
sys.path.append(igrec_dir + "/src/ig_tools/python_utils")
sys.path.append(igrec_dir + "/src/python_pipeline/")
import support
sys.path.append(igrec_dir + "/src/extra/ash_python_utils/")
from ash_python_utils import idFormatByFileName, smart_open, mkdir_p, fastx2fastx, FakeLog
from ig_compress_equal_clusters import parse_cluster_mult


path_to_ig_simulator = igrec_dir + "/../ig_simulator/"
path_to_mixcr = igrec_dir + "/src/extra/aimquast/mixcr/"
path_to_igrec = igrec_dir


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


def multiplex_repertoire(input_file, output_file):
    output_format = idFormatByFileName(output_file)
    input_format = idFormatByFileName(input_file)
    assert output_format == "fasta" or input_format == "fastq"

    with smart_open(input_file) as fh, smart_open(output_file, "w") as fout:
        for record in SeqIO.parse(fh, input_format):
            cluster, mult = parse_cluster_mult(str(record.description))
            for i in xrange(mult):
                record.id = record.description = "antibody_%s_multiplicity_%d_copy_%d" % (cluster, mult, i)
                SeqIO.write(record, fout, output_format)


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


def jit_fx_file(input_file, output_file, error_rate=2, random_errors=True,
                min_error=0, erroneous_site_len=10005000, seed=None):
    import numpy as np
    from Bio import Seq
    import random

    output_format = idFormatByFileName(output_file)
    input_format = idFormatByFileName(input_file)
    random.seed(seed)

    with smart_open(input_file) as fh, smart_open(output_file, "w") as fout:
        for record in SeqIO.parse(fh, input_format):
            n_errors = np.random.poisson(error_rate, 1)[0] if random_errors else error_rate
            if n_errors < min_error:
                n_errors = min_error

            positions = random.sample(range(min(len(record.seq), erroneous_site_len)), n_errors)
            s = list(str(record.seq))
            for pos in positions:
                s[pos] = RC(s[pos])

            if input_format == "fastq":
                phred_quality = record.letter_annotations["phred_quality"]
                record.letter_annotations = {}

            record.seq = Seq.Seq("".join(s))

            if output_format == "fastq":
                if input_format == "fastq":
                    record.letter_annotations["phred_quality"] = phred_quality
                else:
                    record.letter_annotations["phred_quality"] = [random.randint(30, 50) for _ in xrange(len(record))]  # TODO Check it out

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
                    record.letter_annotations["phred_quality"] = [random.randint(30, 50) for _ in xrange(len(record))]  # TODO Check it out
                    # record.letter_annotations["phred_quality"] = [50] * len(record)

                SeqIO.write(record, fout, output_format)


def simulate_data(input_file, output_dir, log=None,
                  **kwargs):
    import tempfile
    import shutil

    if log is None:
        log = FakeLog()

    mkdir_p(output_dir)

    temp_dir = tempfile.mkdtemp()
    run_igrec(input_file, temp_dir, remove_tmp=False, tau=1)  # Run IgReC for VJF output

    input_file = temp_dir + "/vj_finder/cleaned_reads.fa"

    simulated_repertoire_to_rcm(input_file, "%s/ideal_final_repertoire.rcm" % output_dir)

    simulated_repertoire_to_final_repertoire(input_file, "%s/ideal_final_repertoire.fa" % output_dir)

    args = {"path": igrec_dir,
            "repertoire": output_dir + "/ideal_final_repertoire.fa",
            "rcm": output_dir + "/ideal_final_repertoire.rcm"}
    support.sys_call("%(path)s/py/ig_compress_equal_clusters.py %(repertoire)s %(repertoire)s -r %(rcm)s" % args,
                     log=log)

    # TODO factor this stage
    jit_fx_file(input_file, "%s/merged_reads.fa" % output_dir, **kwargs)

    shutil.rmtree(temp_dir)


def run_ig_simulator(output_dir, log=None,
                     chain="HC", num_bases=100, num_mutated=1000, reprtoire_size=5000):
    if log is None:
        log = FakeLog()

    assert chain in ["HC", "LC"]

    args = {"path": path_to_ig_simulator,
            "output_dir": output_dir,
            "chain": chain,
            "num_bases": num_bases,
            "num_mutated": num_mutated,
            "reprtoire_size": reprtoire_size}

    support.sys_call("%(path)s/ig_simulator.py --chain-type %(chain)s --num-bases %(num_bases)d --num-mutated %(num_mutated)d --repertoire-size %(reprtoire_size)d -o %(output_dir)s --skip-drawing" % args,
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


def run_igrec(input_file, output_dir, log=None,
              tau=4,
              min_fillin=0.6,
              loci="all", threads=16, additional_args="",
              remove_tmp=True):
    if log is None:
        log = FakeLog()

    args = {"path": path_to_igrec,
            "tau": tau,
            "min_fillin": min_fillin,
            "loci": loci,
            "threads": threads,
            "input_file": input_file,
            "output_dir": output_dir,
            "additional_args": additional_args}
    support.sys_call("%(path)s/igrec.py --tau=%(tau)d --min-fillin=%(min_fillin)f -t %(threads)d --loci %(loci)s -s %(input_file)s -o %(output_dir)s %(additional_args)s" % args,
                     log=log)
    if remove_tmp:
        import os.path
        if os.path.isfile(output_dir + "/vj_finder"):
            import shutil
            shutil.rmtree(output_dir + "/vj_finder")


def run_mixcr(input_file, output_dir,
              log=None,
              loci="all", enforce_fastq=False,
              threads=16,
              remove_tmp=True):
    if log is None:
        log = FakeLog()

    mkdir_p(output_dir)

    if enforce_fastq and idFormatByFileName(input_file) == "fasta":
        input_file_fq = "%s/input_reads.fq" % output_dir
        fastx2fastx(input_file, input_file_fq)
        input_file = input_file_tmp = input_file_fq
    elif idFormatByFileName(input_file) == "fasta":
        input_file_fasta = "%s/input_reads.fasta" % output_dir
        fastx2fastx(input_file, input_file_fasta)
        input_file = input_file_tmp = input_file_fasta
    else:
        input_file_tmp = None

    args = {"path": path_to_mixcr,
            "compress_eq_clusters_cmd": path_to_igrec + "/py/ig_compress_equal_clusters.py",
            "mixcr_cmd": "java -jar %s/mixcr.jar" % path_to_mixcr,
            "threads": threads,
            "input_file": input_file,
            "output_dir": output_dir,
            "loci": loci}

    support.sys_call("%(mixcr_cmd)s align -t %(threads)d -f -g -r %(output_dir)s/align_report.txt --loci %(loci)s --noMerge -OvParameters.geneFeatureToAlign=VTranscript %(input_file)s %(output_dir)s/mixcr.vdjca" % args,
                     log=log)
    support.sys_call("%(mixcr_cmd)s assemble -t %(threads)d -f -r %(output_dir)s/assemble_report.txt -OassemblingFeatures=\"{FR1Begin:FR4Begin}\" %(output_dir)s/mixcr.vdjca %(output_dir)s/mixcr.clns" % args,
                     log=log)
    support.sys_call("%(mixcr_cmd)s exportClones -pf %(path)s/preset.pf -f --no-spaces %(output_dir)s/mixcr.clns %(output_dir)s/mixcr.txt" % args,
                     log=log)
    convert_mixcr_output_to_igrec("%(output_dir)s/mixcr.txt" % args, "%(output_dir)s/mixcr_uncompressed.fa" % args)
    support.sys_call("%(compress_eq_clusters_cmd)s %(output_dir)s/mixcr_uncompressed.fa %(output_dir)s/final_repertoire.fa" % args)

    if remove_tmp:
        if input_file_tmp is not None:
            os.remove(input_file_tmp)

        os.remove(output_dir + "/mixcr.clns")
        os.remove(output_dir + "/mixcr.txt")
        os.remove(output_dir + "/mixcr.vdjca")
        os.remove(output_dir + "/mixcr_uncompressed.fa")


def parse_presto_id(id):
    import re
    id = id.strip()

    m = re.match(r".*DUPCOUNT=(\d+)$", id)

    if m:
        g = m.groups()
        return int(g[0])
    else:
        return None

assert parse_presto_id("@187968_merged_read_MISEQ2:53:000000000-A2BMV:1:1108:11153:7341_1:N:0:TTAGGC|DUPCOUNT=666") == 666


def run_presto(input_file, output_dir,
               log=None,
               remove_tmp=True):
    if log is None:
        log = FakeLog()

    mkdir_p(output_dir)

    # gunzip
    input_file_new = "%s/input_reads.fasta" % output_dir
    fastx2fastx(input_file, input_file_new)

    args = {"input_file": input_file_new,
            "output_dir": output_dir}

    support.sys_call("CollapseSeq.py -s %(input_file)s --outdir %(output_dir)s --outname presto" % args,
                     log=log)

    presto_output = output_dir + "/presto_collapse-unique.fasta"
    repertoire_fa = output_dir + "/final_repertoire.fa"
    with smart_open(presto_output) as fin, smart_open(repertoire_fa, "w") as fout:
        for i, record in enumerate(SeqIO.parse(fin, idFormatByFileName(presto_output))):
            id = record.description
            size = parse_presto_id(id)
            record.id = record.description = "cluster___%d___size___%d" % (i, size)
            SeqIO.write(record, fout, "fasta")

    if remove_tmp:
        os.remove(input_file_new)
        os.remove(presto_output)


if __name__ == "__main__":
    ig_simulator_output_dir = "/tmp/ig_simulator"
    output_dir = igrec_dir + "/aimquast_test_dataset"
    run_ig_simulator(ig_simulator_output_dir)
    simulate_data(ig_simulator_output_dir + "/final_repertoire.fasta", output_dir)

    run_igrec(output_dir + "/merged_reads.fa",
              output_dir + "/igrec_good/")

    run_igrec(output_dir + "/merged_reads.fa",
              output_dir + "/igrec_bad/", additional_args="--create-triv-dec")
