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

def visualize_region_lengths(labeling_df, region, region_name, output_fname):
    region_seq = list(labeling_df[region])
    region_len = [len(s) for s in region_seq if len(s) > 1]
    sns.distplot(region_len, kde = False, rug=False)
    plt.xlabel(region_name + ' length', fontsize = 16)
    plt.ylabel('# ' + region_name + 's', fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlim(0, 100)
    pp = PdfPages(output_fname + ".pdf")
    pp.savefig()
    plt.savefig(output_fname + ".png")
    pp.close()
    plt.clf()
    print region_name + " length distribution was written to " + output_fname + ".pdf and .png"

############################################################################

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
    ax.legend(ncol = 4, loc="upper center", frameon=True, fontsize = 16)
    plt.xlabel(region_name + ' position', fontsize = 16)
    plt.ylabel('Nucleotide %', fontsize = 16)
    plt.xticks(x, x_l, fontsize = 14)
    plt.yticks(fontsize = 14)
    pp = PdfPages(output_fname + ".pdf")
    pp.savefig()
    plt.savefig(output_fname + ".png")
    pp.close()
    plt.clf()
    print region_name + " nucleotide distribution was written to " + output_fname + ".pdf and .png"

############################################################################

np.random.seed(1)
cm = plt.get_cmap('gnuplot2')
aa_colors = []
for i in range(21):
    aa_colors.append(cm(float(i + 1) / 21)) #(np.random.random_sample()))
amino_acids = ['A', 'G', 'L', 'R', 'W', 'N', 'V', 'I', 'P', 'F', 'Y', 'C', 'T', 'S', 'M', 'Q', 'K', 'H', 'D', 'E', '*']

def visualize_largest_group_aa_variability(labeling_df, region, region_name, output_fname):
    region_seq = list(labeling_df[region])
    max_group = get_region_largest_group(region_seq)
    if len(max_group) == 0:
        return
    group_len = len(max_group[0])
    if group_len % 3 != 0:
        print "Largest " + region_name + " is not out-of-frame"
        return
    aa_seqs = [Seq(cdr).translate(to_stop=True) for cdr in max_group]
    aa_list = [dict() for i in range(0, group_len / 3)]
    for aa_seq in aa_seqs:
        for i in range(0, len(aa_seq)):
            if aa_seq[i] not in aa_list[i]:
                aa_list[i][aa_seq[i]] = 0
            aa_list[i][aa_seq[i]] += 1
    aa_num = [len(aa) for aa in aa_list]
    aa_large_abun = []
    aa_large_acid = []
    for aa in aa_list:
        aa = sorted(aa.items(), key=operator.itemgetter(1), reverse = True)
        sum = 0
        for i in aa:
            sum += i[1]
        aa_large_abun.append(float(aa[0][1]) / float(sum) * 100)
        aa_large_acid.append(aa[0][0])
    aa_set = set()
    for aa in aa_large_acid:
        aa_set.add(aa)
    for aa in aa_set:
        x_ = []
        abun_ = []
        for i in range(0, len(aa_large_abun)):
            x_.append(i)
            if aa_large_acid[i] == aa:
                abun_.append(aa_large_abun[i])
            else:
                abun_.append(0)
        sns.barplot(x_, abun_, color = aa_colors[amino_acids.index(aa)])
    plt.xticks(range(0, len(aa_large_abun)), aa_large_acid, fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlabel('The most abundant amino acid', fontsize = 16)
    plt.ylabel('% ' + region_name + 's', fontsize = 16)
    pp = PdfPages(output_fname + ".pdf")
    pp.savefig()
    plt.savefig(output_fname + ".png")
    pp.close()
    plt.clf()
    print region_name + " aa variability was written to " + output_fname + ".pdf and .png"

def output_cdr_stats_for_locus(locus_df, locus_name, column_name, region_name, output_dir):
    visualize_region_lengths(locus_df, column_name, locus_name + " " + region_name,
                             os.path.join(output_dir, locus_name + "_" + region_name + "_length"))
    visualize_largest_region_nucls(locus_df, column_name, locus_name + " " + region_name,
                             os.path.join(output_dir, locus_name + "_" + region_name + "_nucls"))
    visualize_largest_group_aa_variability(locus_df, column_name, locus_name + " " + region_name,
                             os.path.join(output_dir, locus_name + "_" + region_name + "_aa"))

def output_cdrs_stats_for_locus(vj_df, locus_name, output_dir):
    locus_df = vj_df.loc[vj_df['Chain_type'] == locus_name]
    num_records = len(vj_df['Read_name'])
    num_locus_records = len(locus_df['Read_name'])
    if float(num_locus_records) / float(num_records) < .05:
        print "Output contains very low number of " + locus_name + " records. Drawing plots was skipped"
        return
    print "CDR statistics for " + locus_name
    output_cdr_stats_for_locus(locus_df, locus_name, "CDR1_nucls", "CDR1", output_dir)
    output_cdr_stats_for_locus(locus_df, locus_name, "CDR2_nucls", "CDR2", output_dir)
    output_cdr_stats_for_locus(locus_df, locus_name, "CDR3_nucls", "CDR3", output_dir)

def main(df_fname, output_dir):
    if not os.path.exists(df_fname):
        print "File containing CDR details " + df_fname + " was not found"
        sys.exit(1)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    print "== Output CDR statistics"
    vj_df = pd.read_table(df_fname, delim_whitespace = True)
    output_cdrs_stats_for_locus(vj_df, "IGH", output_dir)
    output_cdrs_stats_for_locus(vj_df, "IGK", output_dir)
    output_cdrs_stats_for_locus(vj_df, "IGL", output_dir)
    print ""

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "Invalid input parameters"
        print "python visualize_cdr_stats.py cdr_details.txt output_dir"
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])