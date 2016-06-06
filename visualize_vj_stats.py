import os
import sys
import operator

from string import letters
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mplt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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
    mplt.rcParams.update({'font.size': 18})
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

def visualize_region_stats(labeling_df, region, region_name, output_fname):
    region_seq = list(labeling_df[region])
    region_len = [len(s) for s in region_seq if len(s) > 1]
    sns.distplot(region_len, kde = False, rug=False)
    plt.xlabel('CDR3 length')
    plt.ylabel('# CDR3s')
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()
    print region_name + " length distribution was written to " + output_fname

def main():
    if len(sys.argv) != 3:
        print "Invalid input parameters"
        return
    df = pd.read_table(sys.argv[1], delim_whitespace = True)
    output_dir = sys.argv[2]
    visualize_vj_heatmap(df, os.path.join(output_dir, "vj_heatmap.pdf"))
    #visualize_region_stats(df, "CDR1_nucls", "CDR1", os.path.join(output_dir, "cdr1_length.pdf"))
    #visualize_region_stats(df, "CDR2_nucls", "CDR2", os.path.join(output_dir, "cdr2_length.pdf"))
    visualize_region_stats(df, "CDR3_nucls", "CDR3", os.path.join(output_dir, "cdr3_length.pdf"))

if __name__ == "__main__":
    main()