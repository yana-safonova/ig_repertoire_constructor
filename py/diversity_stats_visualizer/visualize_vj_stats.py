import os
import sys
import operator
import warnings

import matplotlib as mplt
mplt.use('Agg')

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import visualize_cdr_stats

class VJMatrix:
    def __init__(self, v_hits, j_hits):
        self.num_records = len(v_hits)
        self.vj_dict = dict()
        self.v_set = set()
        self.j_set = set()
        for i in range(0, len(v_hits)):
            self.v_set.add(v_hits[i])
            self.j_set.add(j_hits[i])
            vj_key = v_hits[i], j_hits[i]
            if vj_key not in self.vj_dict:
                self.vj_dict[vj_key] = 0
            self.vj_dict[vj_key] += 1
        print str(len(self.vj_dict)) + " VJ pairs were extracted. Pairs are presented by " + str(len(self.v_set)) + \
              " V genes & " + str(len(self.j_set)) + " J genes"

    def _ComputeVJHitsAbundances(self):
        v_abundance = dict()
        j_abundance = dict()
        for vj in self.vj_dict:
            if vj[0] not in v_abundance:
                v_abundance[vj[0]] = 0
            v_abundance[vj[0]] += self.vj_dict[vj]
            if vj[1] not in j_abundance:
                j_abundance[vj[1]] = 0
            j_abundance[vj[1]] += self.vj_dict[vj]
        return v_abundance, j_abundance

    def _FilterVJByThreshold(self, v_abun, j_abun, size_threshold):
        sorted_v = sorted(v_abun.items(), key=operator.itemgetter(1), reverse=True)
        sorted_j = sorted(j_abun.items(), key=operator.itemgetter(1), reverse=True)
        #print sorted_v
        #print sorted_j
        min_v = min(len(sorted_v), 24)
        min_j = min(len(sorted_j), 100)
        v_abun_large = [sorted_v[i][0] for i in range(0, min_v)] #if float(sorted_v[i][1]) / float(self.num_records) < .4]
        j_abun_large = [sorted_j[i][0] for i in range(0, min_j)] #if float(sorted_j[i][1]) / float(self.num_records) < .4]
        #print v_abun_large
        #print j_abun_large
        return v_abun_large, j_abun_large

    def CreateTable(self, size_threshold = 0):
        v_abun, j_abun = self._ComputeVJHitsAbundances()
        v_abun_l, j_abun_l = self._FilterVJByThreshold(v_abun, j_abun, size_threshold)
        sorted_v_abun_l = sorted(v_abun_l)
        sorted_j_abun_l = sorted(j_abun_l)
        table = [[0] * len(sorted_v_abun_l) for j in sorted_j_abun_l]
        #print table
        for vj in self.vj_dict:
            if vj[0] in sorted_v_abun_l and vj[1] in sorted_j_abun_l:
                v_index = sorted_v_abun_l.index(vj[0])
                j_index = sorted_j_abun_l.index(vj[1])
                table[j_index][v_index] = self.vj_dict[vj]
        sorted_j_abun_l.reverse()
        return table, sorted_v_abun_l, sorted_j_abun_l

def visualize_vj_heatmap(labeling_df, output_pdf):
    v_hits = list(labeling_df['V_hit'])
    j_hits = list(labeling_df['J_hit'])
    vj_matrix = VJMatrix(v_hits, j_hits)
    table, v, j = vj_matrix.CreateTable(100)
    mplt.rcParams.update({'font.size': 20})
    f, ax = plt.subplots(figsize=(15, 15))
    sns.heatmap(table, cmap = plt.cm.jet, xticklabels = v, yticklabels = j, ax = ax)
    ax.tick_params(labelsize = 16)
    x = [i + 0.0 for i in range(0, len(v))]
    y = [i + .5 for i in range(0, len(j))]
    plt.xticks(x, v, rotation=60, fontsize=16)
    plt.yticks(y, j, rotation='horizontal', fontsize=16)
    pp = PdfPages(output_pdf)
    pp.savefig()
    pp.close()
    plt.clf()
    print "VJ heatmap was written to " + output_pdf

############################################################################

def get_gene_isotype(gene_record):
    #return np.random.randint(3) % 3
    splits = gene_record.id.split('|')
    return splits[2].split(':')[1]

def output_shms_pos(all_shms_pos, colors):
    pos = []
    labels = []
    cols = []
    for isotype in all_shms_pos:
        if len(all_shms_pos[isotype]) > 0:
            pos.append(all_shms_pos[isotype])
            labels.append(str(isotype))
            cols.append(colors[isotype])
    plt.hist(pos, bins= 30, color = cols, alpha = .5, label = labels)
    plt.legend(loc = 'upper center', ncol = len(pos), fontsize = 16)
    plt.xlabel("Relative position of SHM in V gene segment", fontsize = 16)
    plt.ylabel("# SHMs", fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)

def output_num_shms(num_all_shms, colors):
    nums = []
    cols = []
    for isotype in num_all_shms:
        if len(num_all_shms[isotype]) > 0:
            sns.distplot(num_all_shms[isotype], hist = False, label = str(isotype))
            #nums.append(num_shms[isotype])
            #labels.append(str(isotype))
            #cols.append(colors[isotype])
    #plt.hist(nums, bins= 25, color = cols, alpha = .5, label = labels)
    #plt.legend(loc = 'upper center', ncol = len(nums))
    plt.xlabel("#SHM in V gene segment", fontsize = 16)
    plt.ylabel("# sequences", fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlim(0, 150)
    plt.legend(fontsize = 14)

def get_nucl_sequence_without_gaps(sequence):
    nucl_seq = Seq("")
    for s in sequence:
        if s != '-':
            nucl_seq += s
    return nucl_seq

def visualize_v_mutations_stats(v_alignment_fasta, output_fname):
    input_records = list(SeqIO.parse(open(v_alignment_fasta), 'fasta'))
    all_shms_pos = {'IGH': [], 'IGK': [], 'IGL': []}
    num_all_shms = {'IGH': [], 'IGK': [], 'IGL': []}
    colors = {'IGH': 'b', 'IGK': 'g', 'IGL': 'r'}
    for i in range(0, len(input_records) / 2):
        read = input_records[i * 2]
        gene = input_records[i * 2 + 1]
        isotype = get_gene_isotype(gene)
        cur_num_shms = 0
        for j in range(0, len(read.seq)):
            if read.seq[j] != gene.seq[j]:
                all_shms_pos[isotype].append(float(j) / float(len(read.seq)))
                cur_num_shms += 1
        num_all_shms[isotype].append(cur_num_shms)
    plt.figure(1, figsize=(9, 6))
    plt.subplot(211)
    output_shms_pos(all_shms_pos, colors)
    plt.subplot(212)
    output_num_shms(num_all_shms, colors)
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()
    print "Distribution of SHMs in V was written to " + output_fname

############################################################################
def skip_aa_record(prev_read_name, cur_read_name, prev_read_pos, cur_read_pos):
    return prev_read_name == cur_read_name and prev_read_pos / 3 == cur_read_pos / 3

def aa_pair_is_valid(aa1, aa2):
    return aa1 != '*' and aa1 != '-' and aa2 != '*' and aa2 != '-'

def get_aa_list():
    return ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def visualize_aa_substitutions(shm_df, output_fname):
    read_aa = shm_df['Read_aa']
    gene_aa = shm_df['Gene_aa']
    read_pos = shm_df['Read_pos']
    read_name = shm_df['Read_name']
    dict_aa = dict()
    prev_read_name = ""
    prev_read_pos = -100
    for i in range(0, len(read_aa)):
        if not skip_aa_record(prev_read_name, read_name[i], prev_read_pos, read_pos[i]) \
                and aa_pair_is_valid(read_aa[i], gene_aa[i]):
            aa_pair = gene_aa[i] + read_aa[i]
            if aa_pair not in dict_aa:
                dict_aa[aa_pair] = 0
            dict_aa[aa_pair] += 1
        prev_read_pos = read_pos[i]
        prev_read_name = read_name[i]
    aa_freq = []
    for i in range(0, 20):
        aa_freq.append([0] * 20)
    aa_list = get_aa_list()
    for aa_pair in dict_aa:
        aa_freq[aa_list.index(aa_pair[1])][aa_list.index(aa_pair[0])] = dict_aa[aa_pair]
    fig, ax = plt.subplots()
    sns.heatmap(aa_freq, cmap = plt.cm.jet, xticklabels = aa_list, yticklabels = aa_list, square = True, ax = ax)
    ax.tick_params(labelsize = 14)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12, rotation='horizontal')
    plt.xlabel("From", fontsize = 14)
    plt.ylabel("To", fontsize = 14, rotation='horizontal')
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()
    print "Amino acid substitution heatmap was written to " + output_fname

def visualize_special_shm_positions(shm_df, output_fname):
    read_pos = shm_df['Read_pos']
    read_length = shm_df['Read_length']
    synonymous = shm_df['Is_synonymous']
    stop_codon = shm_df['Is_stop_codon']
    synonymous_pos = []
    stop_codon_pos = []
    for i in range(0, len(read_pos)):
        if synonymous[i] == 1:
            synonymous_pos.append(float(read_pos[i]) / float(read_length[i]))
        elif stop_codon[i] == 1:
            stop_codon_pos.append(float(read_pos[i]) / float(read_length[i]))
    plt.hist([synonymous_pos, stop_codon_pos], bins= 60, color = ['g', 'r'], alpha = .75, label = ['Synonymous SHMs', 'Stop codon SHMs'])
    plt.legend(loc = 'upper center', ncol = 2, fontsize = 16)
    plt.xlabel("Relative position on read", fontsize = 16)
    plt.ylabel("# SHMs", fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()
    print "Distribution of synonymous and stop codon SHMs was written to " + output_fname

############################################################################

def checkout_output_dir(output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

def main(argv):
    warnings.filterwarnings('ignore')
    if len(argv) != 5:
        print "Invalid input parameters"
        print "python visualize_vj_stats.py cdr_details.txt v_alignment.fasta shm_details.txt output_dir"
        return
    vj_df = pd.read_table(argv[1], delim_whitespace = True)
    v_alignment_fa = argv[2]
    output_dir = argv[4]
    checkout_output_dir(output_dir)
    visualize_vj_heatmap(vj_df, os.path.join(output_dir, "vj_heatmap.pdf"))
    visualize_cdr_stats.main(argv[1], output_dir)
    visualize_v_mutations_stats(v_alignment_fa, os.path.join(output_dir, "v_mutations_distribution.pdf"))

    #shms_df = pd.read_table(argv[3], delim_whitespace = True)
    #visualize_aa_substitutions(shms_df, os.path.join(output_dir, "aa_substitution_heatmap.pdf"))
    #visualize_special_shm_positions(shms_df, os.path.join(output_dir, "synonymous_stop_codon_shms.pdf"))

if __name__ == "__main__":
    main(sys.argv)
