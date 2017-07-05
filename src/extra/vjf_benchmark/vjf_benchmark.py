#!/usr/bin/env python2

import os.path
import sys
import csv
from Bio import SeqIO
import time
from argparse import ArgumentParser
import tempfile


current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = current_dir + "/../../../"
sys.path.append(igrec_dir + "/py/utils")
sys.path.append(igrec_dir + "/py/pipeline/")
sys.path.append(igrec_dir + "/py")
from igblast_utils import ParseIgBlastOutput
import support
sys.path.append(igrec_dir + "/py/")
from ash_python_utils import CreateLogger, AttachFileLogger, linear_search, idFormatByFileName, smart_open, md5_file, fq2fa

from simulate import run_vjfinder

class HitTableRowVJF:
    def __init__(self, _type, query_id, subject_id, start, end):
        self.type = _type
        self.query_id = query_id
        self.subject_id = subject_id
        self.start = start
        self.end = end


def parse_vjf_output(filename, readfile):
    from collections import defaultdict

    with smart_open(readfile, "rU") as fh:
        parser = SeqIO.parse(fh, idFormatByFileName(readfile))
        descr_to_ind = { str(record.description).replace(" ", "_"): i for i, record in enumerate(parser) }

    result = defaultdict(dict)
    with open(filename) as csv_file:
        reader = csv.reader(csv_file, delimiter="\t")
        headers = reader.next()

# Read_name    Chain_type    V_hit    V_start_pos    V_end_pos    V_score    J_hit    J_start_pos    J_end_pos    J_score
        id_col = linear_search(headers, "Read_name")
        Vstart_col = linear_search(headers, "V_start_pos")
        Vend_col = linear_search(headers, "V_end_pos")
        Vgene_col = linear_search(headers, "V_hit")
        Jgene_col = linear_search(headers, "J_hit")
        Jstart_col = linear_search(headers, "J_start_pos")
        Jend_col = linear_search(headers, "J_end_pos")
        for line in reader:
            desc = line[id_col]

            Vstart = int(line[Vstart_col])
            Vend = int(line[Vend_col])
            Jstart = int(line[Jstart_col])
            Jend = int(line[Jend_col])

            Vgene = line[Vgene_col]
            # Vgene = Vgene[:Vgene.find(" ")]
            Jgene = line[Jgene_col]
            # Jgene = Jgene[:Jgene.find(" ")]

            ind = descr_to_ind[desc]
            result[desc]["V"] = HitTableRowVJF("V", desc, Vgene, Vstart, Vend)
            result[desc]["J"] = HitTableRowVJF("J", desc, Jgene, Jstart, Jend)
            result[ind] = result[desc]

        return result


class Empty:
    pass


def hash_file(fname, len=7):
    import os.path

    hash = md5_file(fname)
    base = os.path.basename(fname)

    return base.split(".")[0] + "_" + hash[:len]


def parse_command_line():
    parser = ArgumentParser(description="Benchmark VJFinder vs IgBLAST")
    parser.add_argument("input",
                        type=str,
                        help="input FASTA/FASTQ file")
    parser.add_argument("--bad-reads", "-b",
                        type=str,
                        default="",
                        help="output FASTQ file for suspicious reads, <empty> for non-producing (default: <empty>)")
    parser.add_argument("--options", "-o",
                        type=str,
                        default="-t 16",
                        help="additional options for VJFinder (default: %(default)s)")
    parser.add_argument("--rerun-igblast", "-G",
                        action="store_true",
                        help="perform IgBLAST")
    parser.add_argument("--do-not-run-vjfinder", "-F",
                        action="store_true",
                        help="perform VJFinder")
    parser.add_argument("--tmp-file", "-T",
                        type=str,
                        default="",
                        help="temporary FASTA file used as IgBLAST input, <empty> for automatically generated name (default: <empty>)")
    parser.add_argument("--storage-dir", "-S",
                        type=str,
                        default=".",
                        help="storage directory for cached IgBlast output (default: %(default)s)")
    parser.add_argument("--log", "-L",
                        type=str,
                        default="",
                        help="log file name, <empty> for non-producing (default: <empty>)")

    args = parser.parse_args()

    if not args.tmp_file:
        args.tmp_file = tempfile.mkstemp(suffix=".fa", prefix="vjf_benchmarking_")[1]

    args.workdir = os.path.dirname(os.path.realpath(__file__))

    args.path = args.workdir + "/../../../"

    # args.germline_J_file = args.path + "/data/germline/human/IG/IGHJ-allP.fa"
    args.germline_J_file = args.path + "/data/germline/human/IG/IGHJ.fa"

    return args


def get_vjf_output(args):
    args.vjfinder_output = args.storage_dir + "/" + args.input_hash + "_vjf"
    if not args.do_not_run_vjfinder or not os.path.exists(args.vjfinder_output):
        vjf_time = time.time()
        # support.sys_call("%(path)s/build/release/bin/ig_kplus_vj_finder --db-directory=%(path)s/data/germline -i %(input)s -o %(vjfinder_output)s --separator=tab --loci=IGH --organism=human %(options)s" % args.__dict__, log)
        run_vjfinder(args.input, args.vjfinder_output, loci="IGH", additional_args="  --pseudogenes=false",
                     log=log)
        vjf_time = time.time() - vjf_time
        log.info("VJFinder time: %fs" % vjf_time)

    return parse_vjf_output("%s/alignment_info.csv" % args.vjfinder_output, args.input)


def get_igblast_output(args):
    args.input_hash = hash_file(args.input)
    args.igblast_output = args.storage_dir + "/" + args.input_hash + ".blast"

    if args.rerun_igblast or not os.path.exists(args.igblast_output + ".gz"):
        log.info("IgBLAST output will be written to " + args.igblast_output + ".gz")
        fq2fa(args.input, args.tmp_file)
        igblast_time = time.time()
        support.sys_call("bash %(workdir)s/blast.sh %(tmp_file)s %(igblast_output)s 2> /dev/null" % args.__dict__, log)
        igblast_time = time.time() - igblast_time
        os.unlink(args.tmp_file)
        support.sys_call("gzip %s --force" % args.igblast_output, log)
        log.info("IgBLAST time: %fs" % igblast_time)

    blast = ParseIgBlastOutput(args.igblast_output + ".gz", log, smart_open)
    # Normalize blast_blocks
    return [line.hit_table for line in blast.blocks]


def benchmark_stats(igblast_hits, vjf_hits, germline_J_map):
    align_info_list = []
    for i, line in enumerate(igblast_hits):
        genes = {}
        for row in line:
            genes[row.type] = row
        line.genes = genes

        align_info = Empty()

        # if "V" not in genes or "J" not in genes or genes["V"].evalue > 0.001 or genes["J"].evalue > 100005000:
        if ("V" not in genes) or (genes["V"].evalue > 0.001):
            align_info.contamination = True
            align_info.evalue = genes["V"].evalue if "V" in genes else 999999
        else:
            align_info.contamination = False
            align_info.evalue = genes["V"].evalue

        align_info.vjf_identified = (i in vjf_hits)
        if align_info.vjf_identified:
            vhit = vjf_hits[i]["V"]
            jhit = vjf_hits[i]["J"]
            align_info.vjf_start = vhit.start
            align_info.vjf_end = jhit.end

        if "V" not in genes:
            align_info.igb_start = -999999
        if "J" not in genes:
            align_info.igb_end = -999999
        if "V" in genes:
            align_info.igb_start = genes["V"].q_start - genes["V"].s_start + 1
        if "J" in genes:
            align_info.igb_end = genes["J"].q_end + len(germline_J_map[genes["J"].subject_id]) - genes["J"].s_end

        align_info_list.append(align_info)

    stats = Empty()
    stats.contaminations = 0
    stats.identified_contaminations = 0
    stats.discarded_ig_reads = 0
    stats.vok, stats.jok, stats.vjok = 0, 0, 0
    stats.all = len(align_info_list)
    stats.aligned = 0
    stats.missed_conts = []
    stats.bad_reads = []

    for i, align_info in enumerate(align_info_list):
        if align_info.contamination:
            stats.contaminations += 1
        if align_info.contamination and not align_info.vjf_identified:
            stats.identified_contaminations += 1
        if not align_info.contamination and not align_info.vjf_identified:
            stats.discarded_ig_reads += 1
        if align_info.contamination and align_info.vjf_identified:
            stats.missed_conts.append(i)

        if align_info.vjf_identified and not align_info.contamination:
            stats.aligned += 1
            if align_info.vjf_start == align_info.igb_start:
                stats.vok += 1
            if align_info.vjf_end == align_info.igb_end:
                stats.jok += 1
            if align_info.vjf_start == align_info.igb_start and align_info.vjf_end == align_info.igb_end:
                stats.vjok += 1
            if align_info.vjf_start != align_info.igb_start or align_info.vjf_end != align_info.igb_end:
                stats.bad_reads.append(i)

    return stats


if __name__ == "__main__":
    args = parse_command_line()

    log = CreateLogger("VJF benchmark")
    if args.log:
        AttachFileLogger(log, args.log)

    with open(args.germline_J_file, "rU") as fh:
        germline_J_parser = SeqIO.parse(fh, "fasta")
        germline_J_map = { str(record.id): str(record.seq) for record in germline_J_parser }

    igblast_hits = get_igblast_output(args)
    vjf_hits = get_vjf_output(args)

    with smart_open(args.input, "rU") as fh:
        parser = SeqIO.parse(fh, idFormatByFileName(args.input))
        reads = list(parser)

    assert len(reads) == len(igblast_hits)
    ids = [str(record.description) for record in reads]
    assert len(set(ids)) == len(reads)

    stats = benchmark_stats(igblast_hits, vjf_hits, germline_J_map)

    if args.bad_reads:
        with smart_open(args.bad_reads, "w") as f:
            SeqIO.write([reads[_] for _ in stats.bad_reads], f, idFormatByFileName(args.bad_reads))
        log.info("Bad reads were written to " + args.bad_reads)

    log.info("Overall reads %d" % stats.all)
    log.info("Identified contaminations %d from %d" % (stats.identified_contaminations, stats.contaminations))
    log.info("Missed contaminations count %d " % len(stats.missed_conts))
    log.info("Discarded Ig reads %d " % stats.discarded_ig_reads)
    log.info("Aligned %d" % stats.aligned)
    log.info("V OK, J OK, VJ OK %d %d %d" % (stats.vok, stats.jok, stats.vjok))
    log.info("OK rate %f" % (float(stats.vjok) / float(stats.aligned)))
    log.info("Ill rate %f" % (1. - float(stats.vjok) / float(stats.aligned)))
    log.info("Missed contaminations: %s" % str(stats.missed_conts))
