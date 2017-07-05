#!/usr/bin/env python2

import os
import sys
from Bio import SeqIO

current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = current_dir + "/../"
sys.path.append(igrec_dir + "/py/utils")
sys.path.append(igrec_dir + "/py/pipeline/")
import support
sys.path.append(igrec_dir + "/py/")
from ash_python_utils import idFormatByFileName, smart_open, mkdir_p, fastx2fastx, FakeLog
from ig_compress_equal_clusters import parse_cluster_mult


path_to_ig_simulator = igrec_dir + "/../ig_simulator/"
path_to_ig_simulator_tcr = igrec_dir + "/../ig_simulator_tcr/"
path_to_mixcr = igrec_dir + "/src/extra/tools/mixcr-1.7"
path_to_mixcr2 = igrec_dir + "/src/extra/tools/mixcr-2.0"
path_to_igrec = igrec_dir
path_to_igrec_old = igrec_dir + "/../ig_repertoire_constructor_old"


class Timer:
    def __init__(self):
        import time
        self.start = time.time()

    def delta(self):
        import time
        return time.time() - self.start

    def stamp(self, filename):
        delta = self.delta()
        with smart_open(filename, "w") as f:
            f.write("%f\n" % delta)

        return delta


def parse_final_repertoire_id(id):
    import re
    id = id.strip()

    m = re.match(r"^antibody_([a-zA-Z_0-9]+)_multiplicity_(\d+)_copy_(\d+)$", id)

    if m:
        g = m.groups()
        return g[0].strip(), int(g[1]), int(g[2])
    else:
        return None


assert parse_final_repertoire_id("antibody_1_multiplicity_1_copy_1") == ("1", 1, 1)
assert parse_final_repertoire_id("antibody_NNAAT_multiplicity_1_copy_1") == ("NNAAT", 1, 1)


def parse_abvitro_assembled_header(id):
    import re
    id = id.strip()

    m = re.match(r"^([ACGTN]+)\|CONSCOUNT=(\d+),.*$", id)

    if m:
        g = m.groups()
        return g[0].strip(), int(g[1])
    else:
        return None

assert parse_abvitro_assembled_header("NNATCACTTATAATCCT|CONSCOUNT=777,1|PRCONS=p5-hIGLC_bs-0") == ("NNATCACTTATAATCCT", 777)


def convert_abvitro_to_repertoire(input_file, output_file):
    output_format = idFormatByFileName(output_file)
    input_format = idFormatByFileName(input_file)
    assert output_format == "fasta" or input_format == "fastq"

    with smart_open(input_file) as fh, smart_open(output_file, "w") as fout:
        for record in SeqIO.parse(fh, input_format):
            cluster, mult = parse_abvitro_assembled_header(str(record.description))
            record.id = record.description = "cluster___%s___size___%d" % (cluster, mult)
            SeqIO.write(record, fout, output_format)


def multiplex_repertoire(input_file, output_file):
    output_format = idFormatByFileName(output_file)
    input_format = idFormatByFileName(input_file)
    assert output_format == "fasta" or input_format == "fastq"

    with smart_open(input_file) as fh, smart_open(output_file, "w") as fout:
        for record in SeqIO.parse(fh, input_format):
            cluster, mult = parse_cluster_mult(str(record.description))
            for i in xrange(1, mult + 1):
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
    print seed
    random.seed(seed)
    np.random.seed(seed)

    print np.random.ranf(1)

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
                record.id = record.description = "cluster___%s___size___%d" % (cluster, size)
                record.letter_annotations = {}

                if output_format == "fastq":
                    record.letter_annotations["phred_quality"] = [random.randint(30, 50) for _ in xrange(len(record))]  # TODO Check it out
                    # record.letter_annotations["phred_quality"] = [50] * len(record)

                SeqIO.write(record, fout, output_format)


def simulate_data_wo_errors(input_file, output_dir, log=None):
    import tempfile
    import shutil

    if log is None:
        log = FakeLog()

    mkdir_p(output_dir)

    temp_dir = tempfile.mkdtemp()
    run_igrec(input_file, temp_dir, remove_tmp=False, tau=1)  # Run IgReC for VJF output

    input_file = temp_dir + "/vj_finder/cleaned_reads.fa"

    simulated_repertoire_to_rcm(input_file, "%s/final_repertoire.rcm" % output_dir)

    simulated_repertoire_to_final_repertoire(input_file, "%s/final_repertoire.fa.gz" % output_dir)

    args = {"path": igrec_dir,
            "repertoire": output_dir + "/final_repertoire.fa.gz",
            "rcm": output_dir + "/final_repertoire.rcm"}
    support.sys_call("%(path)s/py/ig_compress_equal_clusters.py %(repertoire)s %(repertoire)s -r %(rcm)s" % args,
                     log=log)

    fastx2fastx(input_file,
                output_dir + "/error_free_reads.fa.gz")

    shutil.rmtree(temp_dir)


def simulate_data_from_dir(input_dir, output_dir, log=None,
                           **kwargs):
    if log is None:
        log = FakeLog()

    mkdir_p(output_dir)

    jit_fx_file(input_dir + "/error_free_reads.fa.gz",
                "%s/input_reads.fa.gz" % output_dir, **kwargs)


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

    simulated_repertoire_to_rcm(input_file, "%s/final_repertoire.rcm" % output_dir)

    simulated_repertoire_to_final_repertoire(input_file, "%s/final_repertoire.fa.gz" % output_dir)

    args = {"path": igrec_dir,
            "repertoire": output_dir + "/final_repertoire.fa.gz",
            "rcm": output_dir + "/final_repertoire.rcm"}
    support.sys_call("%(path)s/py/ig_compress_equal_clusters.py %(repertoire)s %(repertoire)s -r %(rcm)s" % args,
                     log=log)

    # TODO factor this stage
    jit_fx_file(input_file, "%s/input_reads.fa.gz" % output_dir, **kwargs)

    shutil.rmtree(temp_dir)


def run_ig_simulator(output_dir, log=None,
                     chain="HC", num_bases=100, num_mutated=1000, repertoire_size=5000,
                     tcr=False):
    if log is None:
        log = FakeLog()

    assert chain in ["HC", "LC"]

    args = {"path": path_to_ig_simulator if not tcr else path_to_ig_simulator_tcr,
            "output_dir": output_dir,
            "chain": chain,
            "num_bases": num_bases,
            "num_mutated": num_mutated,
            "repertoire_size": repertoire_size}

    timer = Timer()
    cmd = "%(path)s/ig_simulator.py --chain-type %(chain)s --num-bases %(num_bases)d --num-mutated %(num_mutated)d --repertoire-size %(repertoire_size)d -o %(output_dir)s --skip-drawing" % args
    if tcr:
        vgenes = igrec_dir + "/data/germline/human/TCR/TRBV.fa"
        jgenes = igrec_dir + "/data/germline/human/TCR/TRBJ.fa"
        dgenes = igrec_dir + "/data/germline/human/TCR/TRBD.fa"
        cmd += " --vgenes=" + vgenes + " --jgenes=" + jgenes + " --dgenes=" + dgenes
    support.sys_call(cmd,
                     log=log)
    timer.stamp(output_dir + "/time.txt")


def convert_mixcr_output_to_igrec(input_file, output_file):
    with smart_open(input_file) as fh, smart_open(output_file, "w") as fout:
        # Skip header
        fh.next()

        for i, line in enumerate(fh):
            seq, size = line.strip().split()
            size = int(size)
            fout.write(">cluster___%d___size___%d\n" % (i, size))
            fout.write(seq + "\n")


def convert_mixcr2_output_to_igrec(input_file, output_file, initial_reads, output_rcm):
    with smart_open(initial_reads) as fh:
        record_ids = [str(record.description) for record in SeqIO.parse(fh, idFormatByFileName(initial_reads))]

    targets = [None] * len(record_ids)
    with smart_open(input_file) as fh, smart_open(output_file, "w") as fout:
        # Skip header
        fh.next()

        for i, line in enumerate(fh):
            seq, size, ids = line.strip().split("\t")
            ids = ids.strip().split(",")
            ids = map(int, ids)
            for id in ids:
                targets[id] = i
            size = int(size)
            assert size <= len(ids)  # WHY?????????????
            # if size != len(ids):
            #     print size
            #     print ids
            size = len(ids)
            fout.write(">cluster___%d___size___%d\n" % (i, size))
            fout.write(seq + "\n")

        empty_num = max(target for target in targets if target is not None) + 1
        # print empty_num
        with smart_open(initial_reads) as fh:
            for j, record in enumerate(SeqIO.parse(fh, idFormatByFileName(initial_reads))):
                if targets[j] is None:
                    targets[j] = empty_num
                    empty_num += 1
                    fout.write(">cluster___%d___size___%d\n" % (targets[j], 1))
                    fout.write(str(record.seq) + "\n")

    with smart_open(output_rcm, "w") as rcm:
        for id, target_cluster in zip(record_ids, targets):
            assert target_cluster is not None
            rcm.write("%s\t%d\n" % (id, target_cluster))


def dict2class(d):
    class Empty:
        pass

    res = Empty()
    res.__dict__ = d
    return res


def run_vjfinder(input_file, output_dir, log=None,
                 loci="all", threads=16, additional_args="",
                 remove_tmp=False):
    if log is None:
        log = FakeLog()

    import os.path
    import os

    args = {"path": path_to_igrec,
            "loci": loci,
            "threads": threads,
            "input_file": input_file,
            "output_dir": output_dir,
            "organism": "human",
            "path_to_germline": igrec_dir + "/data/germline",
            "additional_args": additional_args}
    args = dict2class(args)

    command_line = args.path + "/build/release/bin/vj_finder" + \
        " -i " + os.path.abspath(args.input_file) + \
        " -o " + os.path.abspath(args.output_dir) + \
        " --db-directory " + os.path.abspath(args.path_to_germline) + \
        " -t " + str(args.threads) + \
        " --loci " + args.loci + \
        " --organism " + args.organism + " " + args.additional_args
    cwd = os.getcwd()
    os.chdir(igrec_dir)
    timer = Timer()
    support.sys_call(command_line, log=log)
    timer.stamp(output_dir + "/time.txt")
    os.chdir(cwd)
    if remove_tmp:
        import os.path
        if os.path.isfile(output_dir):
            import shutil
            shutil.rmtree(output_dir)


def rmdir(dir):
    import os.path
    if os.path.isdir(dir):
        import shutil
        shutil.rmtree(dir)


def run_igrec_old(input_file, output_dir, log=None,
                  tau=3,
                  threads=16, additional_args="",
                  remove_tmp=True):
    if log is None:
        log = FakeLog()

    output_dir = os.path.abspath(output_dir)
    input_file = os.path.abspath(input_file)
    args = {"path": path_to_igrec_old,
            "tau": tau,
            "threads": threads,
            "input_file": input_file,
            "output_dir": output_dir,
            "additional_args": additional_args}
    timer = Timer()
    cwd = os.getcwd()
    os.chdir(path_to_igrec_old)
    support.sys_call("%(path)s/ig_repertoire_constructor.py --tau=%(tau)d -t %(threads)d -s %(input_file)s -o %(output_dir)s %(additional_args)s" % args,
                     log=log)
    os.chdir(cwd)
    timer.stamp(output_dir + "/time.txt")

    # Rename output
    os.rename(output_dir + "/constructed_repertoire.clusters.fa",
              output_dir + "/final_repertoire.fa")
    os.rename(output_dir + "/constructed_repertoire.rcm",
              output_dir + "/final_repertoire.rcm")

    if remove_tmp:
        rmdir(output_dir + "/configs")
        rmdir(output_dir + "/saves")
        rmdir(output_dir + "/temp_files")
        rmdir(output_dir + "/hamming_graphs_tau_%d" % tau)


def run_igrec(input_file, output_dir, log=None,
              tau=4,
              min_fillin=0.6,
              loci="all", threads=16, additional_args="",
              min_sread_size=5,
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
            "min_sread_size": min_sread_size,
            "additional_args": additional_args}
    timer = Timer()
    support.sys_call("%(path)s/igrec.py --tau=%(tau)d --min-fillin=%(min_fillin)f -t %(threads)d --loci %(loci)s -s %(input_file)s -o %(output_dir)s --min-sread-size %(min_sread_size)d %(additional_args)s" % args,
                     log=log)
    timer.stamp(output_dir + "/time.txt")
    if remove_tmp:
        rmdir(output_dir + "/vj_finder")


def run_mixcr(input_file, output_dir,
              log=None,
              loci="all", enforce_fastq=False,
              threads=16,
              remove_tmp=True,
              species="hsa",
              version=1):
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

    path = path_to_mixcr if version == 1 else path_to_mixcr2
    args = {"path": path,
            "compress_eq_clusters_cmd": path_to_igrec + "/py/ig_compress_equal_clusters.py",
            "mixcr_cmd": "java -jar %s/mixcr.jar" % path,
            "threads": threads,
            "input_file": input_file,
            "output_dir": output_dir,
            "species": species,
            "loci": loci,
            "loci_arg": "loci" if version == 1 else "chains"}

    timer = Timer()
    support.sys_call("%(mixcr_cmd)s align -t %(threads)d -f -g -r %(output_dir)s/align_report.txt --%(loci_arg)s %(loci)s --noMerge -OvParameters.geneFeatureToAlign=VTranscript --species %(species)s %(input_file)s %(output_dir)s/mixcr.vdjca" % args,
                     log=log)
    support.sys_call("%(mixcr_cmd)s assemble -t %(threads)d -f -r %(output_dir)s/assemble_report.txt -OassemblingFeatures=\"{FR1Begin:FR4Begin}\" %(output_dir)s/mixcr.vdjca %(output_dir)s/mixcr.clns" % args,
                     log=log)
    support.sys_call("%(mixcr_cmd)s exportClones -sequence -count -f --no-spaces %(output_dir)s/mixcr.clns %(output_dir)s/mixcr.txt" % args,
                     log=log)
    timer.stamp(output_dir + "/time.txt")

    args["features"] = "-count -sequence -nFeature CDR3 -vHit -jHit -vAlignment -jAlignment -aaFeature CDR3"
    support.sys_call("%(mixcr_cmd)s exportClones %(features)s -f --no-spaces %(output_dir)s/mixcr.clns %(output_dir)s/features.txt" % args,
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


def run_mixcr2(input_file, output_dir,
               log=None,
               loci="all", enforce_fastq=False,
               threads=16,
               remove_tmp=True,
               species="hsa",
               region_from="FR1Begin", region_to="FR4Begin"):
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

    path = path_to_mixcr2
    args = {"path": path,
            "compress_eq_clusters_cmd": path_to_igrec + "/py/ig_compress_equal_clusters.py",
            "mixcr_cmd": "java -jar %s/mixcr.jar" % path,
            "threads": threads,
            "input_file": input_file,
            "output_dir": output_dir,
            "species": species,
            "loci": loci,
            "from": region_from,
            "to": region_to,
            "loci_arg": "chains"}

    # support.sys_call("%(mixcr_cmd)s align -t %(threads)d -f -g -r %(output_dir)s/align_report.txt --%(loci_arg)s %(loci)s --noMerge --species %(species)s %(input_file)s %(output_dir)s/mixcr.vdjca" % args,
    timer = Timer()
    #                  log=log)
    support.sys_call("%(mixcr_cmd)s align -p kaligner2 --species %(species)s -t %(threads)d -f -g -r %(output_dir)s/align_report.txt --noMerge --%(loci_arg)s %(loci)s -OreadsLayout=Collinear -OvParameters.geneFeatureToAlign=VTranscript -OallowPartialAlignments=true %(input_file)s %(output_dir)s/mixcr.vdjca" % args,
                     log=log)
    # support.sys_call("%(mixcr_cmd)s assemble -p default_affine -OassemblingFeatures=VDJRegion -OseparateByC=true -OqualityAggregationType=Average -OclusteringFilter.specificMutationProbability=1E-5 -OmaxBadPointsPercent=0 -t %(threads)d -r %(output_dir)s/assemble_report.txt --index %(output_dir)s/index_file %(output_dir)s/mixcr.vdjca %(output_dir)s/mixcr.clns" % args,
    # support.sys_call("%(mixcr_cmd)s assemble -f -p default_affine -OassemblingFeatures=VDJRegion -OseparateByC=true -OqualityAggregationType=Average -OclusteringFilter.specificMutationProbability=1E-5 -OmaxBadPointsPercent=0 -r %(output_dir)s/assemble_report.txt --index %(output_dir)s/index_file %(output_dir)s/mixcr.vdjca %(output_dir)s/mixcr.clns" % args,
    #                  log=log)
    # support.sys_call("%(mixcr_cmd)s assemble -t %(threads)d -f -r %(output_dir)s/assemble_report.txt --index %(output_dir)s/index_file %(output_dir)s/mixcr.vdjca %(output_dir)s/mixcr.clns" % args,
    #                  log=log)
    support.sys_call("%(mixcr_cmd)s assemble -t %(threads)d -f -r %(output_dir)s/assemble_report.txt --index %(output_dir)s/index_file -OassemblingFeatures=\"{%(from)s:%(to)s}\" %(output_dir)s/mixcr.vdjca %(output_dir)s/mixcr.clns" % args,
                     log=log)
    args["small_features"] = "-sequence -count -readIds %(output_dir)s/index_file" % args
    support.sys_call("%(mixcr_cmd)s exportClones %(small_features)s -f --no-spaces %(output_dir)s/mixcr.clns %(output_dir)s/mixcr.txt" % args,
                     log=log)
    timer.stamp(output_dir + "/time.txt")

    args["features"] = "-count -sequence -nFeature CDR3 -vHit -jHit -vAlignment -jAlignment -aaFeature CDR3 -readIds %(output_dir)s/index_file" % args
    support.sys_call("%(mixcr_cmd)s exportClones %(features)s -f --no-spaces %(output_dir)s/mixcr.clns %(output_dir)s/features.txt" % args,
                     log=log)
    # convert_mixcr_output_to_igrec("%(output_dir)s/mixcr.txt" % args, "%(output_dir)s/mixcr_uncompressed.fa" % args)

    convert_mixcr2_output_to_igrec("%(output_dir)s/mixcr.txt" % args,
                                   "%(output_dir)s/mixcr_uncompressed.fa" % args,
                                   input_file,
                                   "%(output_dir)s/mixcr_uncompressed.rcm" % args)
    support.sys_call("%(compress_eq_clusters_cmd)s %(output_dir)s/mixcr_uncompressed.fa %(output_dir)s/final_repertoire.fa -r %(output_dir)s/mixcr_uncompressed.rcm -R %(output_dir)s/final_repertoire.rcm" % args)

    if remove_tmp:
        if input_file_tmp is not None:
            os.remove(input_file_tmp)

        os.remove(output_dir + "/align_report.txt")
        os.remove(output_dir + "/assemble_report.txt")
        os.remove(output_dir + "/mixcr.clns")
        os.remove(output_dir + "/mixcr.txt")
        os.remove(output_dir + "/features.txt")
        os.remove(output_dir + "/mixcr.vdjca")
        os.remove(output_dir + "/mixcr_uncompressed.fa")
        os.remove(output_dir + "/mixcr_uncompressed.rcm")
        os.remove(output_dir + "/index_file")


def run_mixcr2_alignment_only(input_file, output_dir,
                              log=None,
                              loci="all", enforce_fastq=False,
                              threads=16,
                              remove_tmp=True,
                              species="hsa"):
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

    path = path_to_mixcr2
    args = {"path": path,
            "compress_eq_clusters_cmd": path_to_igrec + "/py/ig_compress_equal_clusters.py",
            "mixcr_cmd": "java -jar %s/mixcr.jar" % path,
            "threads": threads,
            "input_file": input_file,
            "output_dir": output_dir,
            "species": species,
            "loci": loci,
            "loci_arg": "chains"}

    # support.sys_call("%(mixcr_cmd)s align -t %(threads)d -f -g -r %(output_dir)s/align_report.txt --%(loci_arg)s %(loci)s --noMerge --species %(species)s %(input_file)s %(output_dir)s/mixcr.vdjca" % args,
    #                  log=log)
    timer = Timer()
    support.sys_call("%(mixcr_cmd)s align -p kaligner2 --species %(species)s -t %(threads)d -f -g -r %(output_dir)s/align_report.txt --noMerge --%(loci_arg)s %(loci)s -OreadsLayout=Collinear -OvParameters.geneFeatureToAlign=VTranscript -OallowPartialAlignments=true %(input_file)s %(output_dir)s/mixcr.vdjca" % args,
                     log=log)
    timer.stamp(output_dir + "/time.txt")

    if remove_tmp:
        if input_file_tmp is not None:
            os.remove(input_file_tmp)

        os.remove(output_dir + "/align_report.txt")
        os.remove(output_dir + "/mixcr.vdjca")


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

    timer = Timer()
    support.sys_call("CollapseSeq.py -s %(input_file)s --outdir %(output_dir)s --outname presto" % args,
                     log=log)
    timer.stamp(output_dir + "/time.txt")

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
    output_dir = igrec_dir + "/test_dataset/igquast"
    run_ig_simulator(ig_simulator_output_dir)
    simulate_data(ig_simulator_output_dir + "/final_repertoire.fasta", output_dir)

    run_igrec(output_dir + "/input_reads.fa",
              output_dir + "/igrec_good/")

    run_igrec(output_dir + "/input_reads.fa",
              output_dir + "/igrec_bad/", additional_args="--create-triv-dec")
