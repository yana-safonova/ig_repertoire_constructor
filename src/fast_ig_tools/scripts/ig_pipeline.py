#!/usr/bin/env python2


def process_file(filename, chain="heavy", threads=16):
    import os.path
    import re

    filename = os.path.abspath(filename)
    dirname = os.path.dirname(filename)
    basename = os.path.basename(filename)

    basename = re.sub(r"\.(gz|bz2)$", "", basename)
    basename = re.sub(r"\.(fa|fq|fasta|fastq)$", "", basename)

    filename_initial = filename
    filename_cropped = "%s/%s.cropped.fa.gz" % (dirname, basename)
    filename_bad = "%s/%s.bad.fa.gz" % (dirname, basename)
    filename_addinfo = "%s/%s.addinfo.csv" % (dirname, basename)
    filename_collapsed = "%s/%s.collapsed.fa.gz" % (dirname, basename)
    filename_graph = "%s/%s.graph" % (dirname, basename)

    os.system("./ig_kplus_vj_finder -i %s -o %s -b %s -a %s --chain=%s -t %d" %
              (filename_initial, filename_cropped, filename_bad, filename_addinfo,
               chain, threads))

    os.system("./ig_trie_compressor -i %s -o %s" %
              (filename_cropped, filename_collapsed))

    os.system("./ig_swgraph_construct -i %s -o %s -t %d" %
              (filename_collapsed, filename_graph, threads))


if __name__ == "__main__":
    process_file("1_SAM13306969_merged.fasta")
    process_file("2_SAM13306970_merged.fasta")
    process_file("3_SAM13306971_merged.fq", "lambda")
    process_file("4_SAM13306972.merged.fastq", "lambda")
    process_file("5_SAM15574990.merged.fastq", "lambda")
    process_file("6_SAM15574989.merged.fastq", "lambda")
    process_file("7_SAM15574987.merged.fastq")
    process_file("8_SAM15574988.merged.fastq")

    process_file("age_ig_s1_R12_new.fastq")
    process_file("age_ig_s2_R12_new.fastq")
    process_file("age_ig_s3_R12_new.fastq")
    process_file("age_ig_s4_R12_new.fastq")
    process_file("age_ig_s5_R12_new.fastq")
    process_file("age_ig_s6_R12_new.fastq")
    process_file("age_ig_s7_R12_new.fastq")
    process_file("age_ig_s8_R12_new.fastq")
    process_file("age_ig_s9_R12_new.fastq")
