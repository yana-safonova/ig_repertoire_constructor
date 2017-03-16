import os
import sys
import operator

import matplotlib as mplt
mplt.use('Agg')

import utils

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from Bio.Seq import Seq

def visualize_region_lengths(labeling_df, region, region_name, output_fname, log):
    region_seq = list(labeling_df[region])
    region_len = [len(s) for s in region_seq if len(s) > 1]
    f, ax = plt.subplots(figsize=(8, 8))
    sns.distplot(region_len, kde = False, rug=False)
    plt.xlabel(region_name + ' length', fontsize = 16)
    plt.ylabel('# ' + region_name + 's', fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlim(0, 100)
    utils.output_figure(output_fname, region_name + " length distribution", log)

############################################################################
def visualize_length_abundance_dist(labeling_df, region, region_name, output_fname, log):
    region_seq = list(labeling_df[region])
    region_dict = dict()
    for seq in region_seq:
        if seq not in region_dict:
            region_dict[seq] = 0
        region_dict[seq] += 1
    abun = [] #np.array()
    lens = [] #np.array()
    for seq in region_dict:
        if region_dict[seq]  == 1:
            continue
        abun.append(region_dict[seq])
        lens.append(len(seq))
    abun = np.asarray(abun)
    lens = np.asarray(lens)
    f, ax = plt.subplots()
    sns.jointplot(abun, lens, size = 6)
    #plt.xlabel(region_name + ' abundance', fontsize = 14)
    #ax.xaxis.set_label_position('top')
    #plt.ylabel(region_name + ' length', fontsize = 14)
    #plt.xticks(fontsize = 14)
    #plt.yticks(fontsize = 14)
    #plt.xlim(-1, abun.max() + 1)
    #plt.ylim(-1, lens.max() + 1)
    utils.output_figure(output_fname, region_name + " joint distribution of abundances & lengths", log)

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

def visualize_largest_region_nucls(labeling_df, region, region_name, output_fname, log):
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
    utils.output_figure(output_fname, region_name + " nucleotide distribution", log)

############################################################################

amino_acids = ['A', 'G', 'L', 'R', 'W', 'N', 'V', 'I', 'P', 'F', 'Y', 'C', 'T', 'S', 'M', 'Q', 'K', 'H', 'D', 'E', '*']

def get_aa_colors():
    aa_colors = []
    for aa in amino_acids:
        hydrophoby = utils.hydrophoby_dict[aa]
        rel_value = float(hydrophoby - utils.hydro_min) / float(utils.hydro_max - utils.hydro_min)
        aa_colors.append(plt.get_cmap('bwr')(rel_value))
    return aa_colors

def visualize_largest_group_aa_variability(labeling_df, region, region_name, output_fname, log):
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
    aa_colors = get_aa_colors()
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
        df = pd.DataFrame({'x': x_, 'y' : abun_})
        sns.barplot(x = 'x', y = 'y', data = df, color = aa_colors[amino_acids.index(aa)])
    plt.xticks(range(0, len(aa_large_abun)), aa_large_acid, fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlabel('The most abundant amino acid', fontsize = 16)
    plt.ylabel('% ' + region_name + 's', fontsize = 16)
    utils.output_figure(output_fname, region_name + " aa variability", log)

def output_cdr_stats_for_locus(locus_df, locus_name, column_name, region_name, output_dir, log):
    visualize_region_lengths(locus_df, column_name, locus_name + " " + region_name,
                             os.path.join(output_dir, locus_name + "_" + region_name + "_length"), log)
    #visualize_length_abundance_dist(locus_df, column_name, locus_name + " " + region_name,
    #                        os.path.join(output_dir, locus_name + "_" + region_name + "_abundance_length"), log)
    visualize_largest_region_nucls(locus_df, column_name, locus_name + " " + region_name,
                             os.path.join(output_dir, locus_name + "_" + region_name + "_nucls"), log)
    visualize_largest_group_aa_variability(locus_df, column_name, locus_name + " " + region_name,
                             os.path.join(output_dir, locus_name + "_" + region_name + "_aa"), log)

def output_cdrs_stats_for_locus(vj_df, locus_name, output_dir, log):
    locus_df = vj_df.loc[vj_df['Chain_type'] == locus_name]
    num_records = len(vj_df['Read_name'])
    num_locus_records = len(locus_df['Read_name'])
    if float(num_locus_records) / float(num_records) < .05 or num_locus_records < 10:
        log.info("Output contains very low number (" + str(num_locus_records) + ") of " + locus_name + " records. Drawing plots was skipped")
        return
    log.info("Visualization of CDR statistics for " + locus_name + " locus")
    output_cdr_stats_for_locus(locus_df, locus_name, "CDR1_nucls", "CDR1", output_dir, log)
    output_cdr_stats_for_locus(locus_df, locus_name, "CDR2_nucls", "CDR2", output_dir, log)
    output_cdr_stats_for_locus(locus_df, locus_name, "CDR3_nucls", "CDR3", output_dir, log)

def main(df_fname, output_dir, log):
    if not os.path.exists(df_fname):
        log.info("File containing CDR details " + df_fname + " was not found")
        sys.exit(1)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    log.info("== Output CDR statistics")
    vj_df = pd.read_table(df_fname, delim_whitespace = True)
    if len(vj_df['Read_name']) == 0:
        log.info("CDR data-frame contains 0 records. CDR visualization will be skipped")
        return
    output_cdrs_stats_for_locus(vj_df, "IGH", output_dir, log)
    output_cdrs_stats_for_locus(vj_df, "IGK", output_dir, log)
    output_cdrs_stats_for_locus(vj_df, "IGL", output_dir, log)
    log.info("")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "Invalid input parameters"
        print "python visualize_cdr_stats.py cdr_details.txt output_dir logger"
        sys.exit(1)
    log = utils.get_logger_by_arg(sys.argv[3], "cdr_visualization")
    main(sys.argv[1], sys.argv[2], log)