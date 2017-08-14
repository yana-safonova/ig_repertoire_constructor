import sys
sys.path.insert(0, "..")

import os
import pandas as pd
from Bio import SeqIO
from spots.spots import hotspots, coldspots

int_kmers = ["TTATA", "GCTCC", "CAACA", "AATAA", "CATCT", "ACAGC", "ATCTA"]
hspots = hotspots()
cspots = coldspots()
for int_kmer in int_kmers:
    assert int_kmer not in hspots
    assert int_kmer not in cspots

kmer_len = 5

if __name__ == '__main__':
    germline_file = "../../../data/germline/human/IG/IGKV.fa"
    with open(germline_file, 'r') as f:
        germline = list(SeqIO.parse(f, "fasta"))
    germline = {str(s.id): str(s.seq).upper() for s in germline}

    annotation_germline_file = "../../../data/annotation/human_v_imgt.txt"
    annotation_germline = pd.read_csv(annotation_germline_file,
                                      sep='\t',
                                      index_col=0, header=None,
                                      names=['gene',
                                             'fr1_start', 'fr1_end', 'cdr1_start', 'cdr1_end',
                                             'fr2_start', 'fr2_end', 'cdr2_start', 'cdr2_end',
                                             'fr3_start', 'fr3_end', 'a', 'b'])

    germline_cdr1, germline_cdr2 = {}, {}
    cdr1_lens, cdr2_lens = [], []
    for gene, seq in germline.iteritems():
        if gene not in annotation_germline.index:
            continue
        cdr1_start = annotation_germline.cdr1_start[gene] - 1 - kmer_len // 2
        cdr1_end = annotation_germline.cdr1_end[gene] + kmer_len // 2
        cdr2_start = annotation_germline.cdr2_start[gene] - 1 - kmer_len // 2
        cdr2_end = annotation_germline.cdr2_end[gene] + kmer_len // 2
        cdr1_start, cdr1_end = int(cdr1_start), int(cdr1_end)
        cdr2_start, cdr2_end = int(cdr2_start), int(cdr2_end)
        cdr1 = seq[cdr1_start:cdr1_end]
        cdr2 = seq[cdr2_start:cdr2_end]
        cdr1_lens.append(len(cdr1))
        cdr2_lens.append(len(cdr2))

        cdr = cdr1 + "XXXXX" + cdr2

        for int_kmer in int_kmers:
            ind = cdr.find(int_kmer)
            ind_xxx = cdr.find("XXXXX")
            if ind == -1:
                continue
            # for i in xrange(max(0, ind - kmer_len // 2), min(len(cdr), ind + kmer_len // 2 + 1)):
            for i in xrange(ind - kmer_len // 2, ind + kmer_len // 2 + 1):
                kmer = cdr[i : i + kmer_len]
                if kmer in hspots:
                    print(int_kmer, kmer, i, "CDR1" if ind < ind_xxx else "CDR2")

    cdr1_lens = list(set(cdr1_lens))
    cdr2_lens = list(set(cdr2_lens))
    print(cdr1_lens)
    print(cdr2_lens)


