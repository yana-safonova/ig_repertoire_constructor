import os
import sys
import operator
import warnings

from string import letters

import matplotlib as mplt
mplt.use('Agg')

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

mplt.use('Agg')

class VJMatrix:
    def __init__(self, v_hits, j_hits):
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
        v_abun_large = [sorted_v[i][0] for i in range(0, min_v)] #[v[0] for v in sorted_v if v[1] >= size_threshold]
        j_abun_large = [sorted_j[i][0] for i in range(0, min_j)] #[j[0] for j in sorted_j if j[1] >= size_threshold]
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
    table, v, j = vj_matrix.CreateTable(2500)
    mplt.rcParams.update({'font.size': 20})
    f, ax = plt.subplots(figsize=(15, 10))
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    #sns.clustermap(table, cmap = plt.cm.coolwarm)
    sns.heatmap(table, cmap = plt.cm.coolwarm, xticklabels = v, yticklabels = j) #plt.cm.Greens)
    x = [i + 0.0 for i in range(0, len(v))]
    y = [i + .5 for i in range(0, len(j))]
    plt.xticks(x, v, rotation=60) #, fontsize=12)
    plt.yticks(y, j, rotation='horizontal') #, fontsize=12)
    pp = PdfPages(output_pdf)
    pp.savefig()
    pp.close()
    plt.clf()
    print "VJ heatmap was written to " + output_pdf

def visualize_region_lengths(labeling_df, region, region_name, output_fname):
    region_seq = list(labeling_df[region])
    region_len = [len(s) for s in region_seq if len(s) > 1]
    sns.distplot(region_len, kde = False, rug=False)
    plt.xlabel(region_name + ' length')
    plt.ylabel('# ' + region_name + 's')
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()
    print region_name + " length distribution was written to " + output_fname

def get_region_largest_group(region_seq):
    len_dict = dict()
    for s in region_seq:
        if len(s) not in len_dict:
            len_dict[len(s)] = list()
        len_dict[len(s)].append(s)
    max_len = 0
    max_group = 0
    for l in len_dict:
        if len(len_dict[l]) > max_group:
            max_group = len(len_dict[l])
            max_len = l
    return len_dict[max_len]

def get_nucls_lists(max_group):
    region_len = len(max_group[0])
    nucl_dict = dict()
    nucl_dict['A'] = [0] * region_len
    nucl_dict['C'] = [0] * region_len
    nucl_dict['G'] = [0] * region_len
    nucl_dict['T'] = [0] * region_len
    for s in max_group:
        for i in range(0, len(s)):
            if s[i] in nucl_dict:
                nucl_dict[s[i]][i] += 1
    for i in range(0, region_len):
        sum = nucl_dict['A'][i] + nucl_dict['C'][i] + nucl_dict['G'][i] + nucl_dict['T'][i]
        nucl_dict['A'][i] = float(nucl_dict['A'][i]) / float(sum) * 100
        nucl_dict['C'][i] = float(nucl_dict['C'][i]) / float(sum) * 100
        nucl_dict['G'][i] = float(nucl_dict['G'][i]) / float(sum) * 100
        nucl_dict['T'][i] = float(nucl_dict['T'][i]) / float(sum) * 100
    nucl_dict['A'] = np.array(nucl_dict['A'])
    nucl_dict['C'] = np.array(nucl_dict['C'])
    nucl_dict['G'] = np.array(nucl_dict['G'])
    nucl_dict['T'] = np.array(nucl_dict['T'])
    return nucl_dict

def visualize_largest_region_nucls(labeling_df, region, region_name, output_fname):
    region_seq = list(labeling_df[region])
    max_group = get_region_largest_group(region_seq)
    if len(max_group) == 0:
        return
    nucl_dict = get_nucls_lists(max_group)
    x = range(0, len(max_group[0]))
    x_l = [str(i) for i in range(1, len(max_group[0]) + 1)]
    acgt = nucl_dict['A'] + nucl_dict['C'] + nucl_dict['G'] + nucl_dict['T']
    cgt = nucl_dict['C'] + nucl_dict['G'] + nucl_dict['T']
    gt = nucl_dict['G'] + nucl_dict['T']
    #sns.set_color_codes("pastel")
    sns.set_color_codes("muted")
    f, ax = plt.subplots(figsize=(15, 6))
    sns.barplot(x=x, y=acgt, label="A", color = 'b')
    sns.barplot(x=x, y=cgt, label="C", color = 'g')
    sns.barplot(x=x, y=gt, label="G", color = 'r')
    sns.barplot(x=x, y=nucl_dict['T'], label="T", color = 'orange')
    ax.legend(ncol = 4, loc="upper center", frameon=True, fontsize = 14)
    plt.xlabel(region_name + ' position', fontsize = 14)
    plt.ylabel('Nucleotide %', fontsize = 14)
    plt.xticks(x, x_l)
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()
    print region_name + " nucleotide distribution was written to " + output_fname

def checkout_output_dir(output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

def main(argv):
    warnings.filterwarnings('ignore')
    if len(argv) != 3:
        print "Invalid input parameters"
        return
    df = pd.read_table(argv[1], delim_whitespace = True)
    output_dir = argv[2]
    checkout_output_dir(output_dir)
    visualize_vj_heatmap(df, os.path.join(output_dir, "vj_heatmap.pdf"))
    visualize_region_lengths(df, "CDR1_nucls", "CDR1", os.path.join(output_dir, "cdr1_length.pdf"))
    visualize_region_lengths(df, "CDR2_nucls", "CDR2", os.path.join(output_dir, "cdr2_length.pdf"))
    visualize_region_lengths(df, "CDR3_nucls", "CDR3", os.path.join(output_dir, "cdr3_length.pdf"))
    visualize_largest_region_nucls(df, "CDR1_nucls", "CDR1", os.path.join(output_dir, "cdr1_nucls.pdf"))
    visualize_largest_region_nucls(df, "CDR2_nucls", "CDR2", os.path.join(output_dir, "cdr2_nucls.pdf"))
    visualize_largest_region_nucls(df, "CDR3_nucls", "CDR3", os.path.join(output_dir, "cdr3_nucls.pdf"))

if __name__ == "__main__":
    main(sys.argv)
