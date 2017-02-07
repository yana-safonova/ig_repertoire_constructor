import os

import pandas as pd
from Bio import SeqIO

from chains.chains import Chains


file_dir = os.path.dirname(os.path.realpath(__file__))
file_dir = os.path.abspath(file_dir)
igrec_path = os.path.normpath(os.path.join(file_dir, "..", "..", ".."))
germline_path = os.path.join(igrec_path, "data", "germline", "human", "IG")
annotation_filename = os.path.join(igrec_path,
                                   "data",
                                   "annotation",
                                   "human_v_imgt.txt")

k = 5

def get_genomic_kmers(chain):
    if chain == Chains.IGH:
        germline_filename = os.path.join(germline_path, "IGHV.fa")
    elif chain == Chains.IGK:
        germline_filename = os.path.join(germline_path, "IGKV.fa")
    elif chain == Chains.IGL:
        germline_filename = os.path.join(germline_path, "IGLV.fa")

    with open(germline_filename, 'r') as f:
        germline = list(SeqIO.parse(f, "fasta"))
    germline = {str(s.id): str(s.seq).upper() for s in germline}

    header = ['gene',
              'fr1_start', 'fr1_end', 'cdr1_start', 'cdr1_end',
              'fr2_start', 'fr2_end', 'cdr2_start', 'cdr2_end',
              'fr3_start', 'fr3_end', 'a', 'b']
    annot = pd.read_csv(annotation_filename,
                        sep='\t',
                        index_col=0, header=None,
                        names=header)
    kmers_fr, kmers_cdr = [], []

    for gene, seq in germline.iteritems():
        if gene not in annot.index:
            continue
        fr1_start  = int(annot.fr1_start [gene] - 1)
        cdr1_start = int(annot.cdr1_start[gene] - 1)
        fr2_start  = int(annot.fr2_start [gene] - 1)
        cdr2_start = int(annot.cdr2_start[gene] - 1)
        fr3_start  = int(annot.fr3_start [gene] - 1)
        fr3_end    = int(annot.fr3_end   [gene] - 1)

        fr1 = seq[fr1_start:cdr1_start]
        fr2 = seq[fr2_start:cdr2_start]
        fr3 = seq[fr3_start:fr3_end]
        cdr1 = seq[cdr1_start:fr2_start]
        cdr2 = seq[cdr2_start:fr3_start]

        kmers_fr += [seq[i : i + k] for i in xrange(len(fr1) - k + 1)]
        kmers_fr += [seq[i : i + k] for i in xrange(len(fr2) - k + 1)]
        kmers_fr += [seq[i : i + k] for i in xrange(len(fr3) - k + 1)]

        kmers_cdr += [seq[i : i + k] for i in xrange(len(cdr1) - k + 1)]
        kmers_cdr += [seq[i : i + k] for i in xrange(len(cdr2) - k + 1)]

    kmers_fr = list(set(kmers_fr))
    kmers_cdr = list(set(kmers_cdr))

    return kmers_fr, kmers_cdr


def get_genomic_kmers_all_chains():
    all_kmers_fr = dict.fromkeys(Chains)
    all_kmers_cdr = dict.fromkeys(Chains)
    for chain in Chains:
        if chain == Chains.IG:
            continue
        kmers_fr, kmers_cdr = get_genomic_kmers(chain)
        all_kmers_fr[chain] = kmers_fr
        all_kmers_cdr[chain] = kmers_cdr
    return all_kmers_fr, all_kmers_cdr
