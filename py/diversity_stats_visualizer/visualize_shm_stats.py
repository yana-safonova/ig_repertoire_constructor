import os
import sys

import utils

import pylab
import matplotlib as mplt
mplt.use('Agg')

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches

class AlignedRead:
    def _parse_line(self, line):
        splits = line.split()
        if len(splits) != 6:
            print "Header line " + line + " is not valid"
            sys.exit(1)
        self.read_name = splits[0][len("Read_name:"):]
        self.read_len = int(splits[1][len("Read_length:"):])
        self.gene_name = splits[2][len("Gene_name:"):]
        self.gene_len = int(splits[3][len("Gene_length:"):])
        self.segment = splits[4][len("Segment:"):]
        self.chain_type = splits[5][len("Chain_type:"):]

    def __init__(self, line):
        self._parse_line(line)

    def __hash__(self):
        return hash(self.read_name) * hash(self.gene_name)

    def __eq__(self, other):
        return self.read_name == other.read_name and self.gene_name == other.gene_name

    def is_variable(self):
        return self.segment == 'V'

class SHM:
    def __init__(self, line):
        splits = line.split()
        if len(splits) != 9:
            print "Invalid SHM line " + line
        self.type = splits[0]
        self.read_pos = int(splits[1])
        self.gene_pos = int(splits[2])
        self.read_nucl = splits[3]
        self.gene_nucl = splits[4]
        self.read_aa = splits[5]
        self.gene_aa = splits[6]
        self.synonymous = False
        if int(splits[7]) == 1:
            self.synonymous = True
        self.to_stop_codon = False
        if int(splits[8]) == 1:
            self.to_stop_codon = True

    def is_deletion(self):
        return self.type == 'D'

    def is_insertion(self):
        return self.type == 'I'

    def is_substitution(self):
        return self.type == 'S'

class SHMs:
    def __init__(self, df_fname):
        self.shm_dict = dict()
        fhandler = open(df_fname, "r")
        lines = fhandler.readlines()
        current_read = ""
        for i in range(1, len(lines)):
            l = lines[i].strip()
            if l[:len("Read_name:")] == "Read_name:":
                current_read = AlignedRead(l)
                self.shm_dict[current_read] = list()
            else:
                self.shm_dict[current_read].append(SHM(l))

    def __len__(self):
        return len(self.shm_dict)

    def __getitem__(self, item):
        return self.shm_dict[item]

    def __iter__(self):
        for it in self.shm_dict:
            yield it


isotype_colors = {'IGH': 'b', 'IGK': 'g', 'IGL': 'r'}

def output_shms_pos(all_shms_pos, colors, output_prefix, log):
    for isotype in all_shms_pos:
        if len(all_shms_pos[isotype]) < 10:
            continue
        plt.hist(all_shms_pos[isotype], bins = 100, color = colors[isotype], alpha = .75)
        plt.xlabel("#SHM in " + isotype + "V gene segment", fontsize = 16)
        plt.ylabel("# sequences", fontsize = 16)
        plt.xticks(fontsize = 14)
        plt.yticks(fontsize = 14)
        plt.xlim(0, .75)
        output_fname = output_prefix + "_" + isotype + "V_pos"
        utils.output_figure(output_fname, "Distribution of SHM relative positions in " + isotype + "V segments", log)

def output_num_shms(num_all_shms, colors, output_prefix, log):
    pos = []
    labels = []
    cols = []
    for isotype in num_all_shms:
        if len(num_all_shms[isotype]) > 0:
            pos.append(num_all_shms[isotype])
            labels.append(str(isotype))
            cols.append(colors[isotype])
    plt.hist(pos, bins = 50, color = cols, alpha = .75, label = labels)
    plt.legend(loc = 'upper center', ncol = len(pos), fontsize = 16)
    plt.xlabel("# of SHM in V gene segment", fontsize = 16)
    plt.ylabel("# SHMs", fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlim(0, 150)
    output_fname = output_prefix + "_shms_number"
    utils.output_figure(output_fname, "Distribution of # SHMs in V segments", log)

def output_shm_stats_for_isotype(num_shms, shm_pos, isotype, output_prefix, log):
    plt.figure(1)
    plt.subplot(211)
    # plot for SHM positions
    plt.hist(shm_pos, color = isotype_colors[isotype], alpha = .75, bins = 50)
    #cdr_color = "#EFBEBE"
    #plt.gca().add_patch(patches.Rectangle((cdr_positions[isotype]['CDR1'][0], 0),
    #                                      cdr_positions[isotype]['CDR1'][1] - cdr_positions[isotype]['CDR1'][0],
    #                                      max(n) + 2, facecolor= cdr_color, lw = 0))
    #plt.gca().add_patch(patches.Rectangle((cdr_positions[isotype]['CDR2'][0], 0),
    #                                      cdr_positions[isotype]['CDR2'][1] - cdr_positions[isotype]['CDR2'][0],
    #                                      max(n) + 2, facecolor= cdr_color, lw = 0))
    #plt.gca().add_patch(patches.Rectangle((cdr_positions[isotype]['CDR3'][0], 0),
    #                                      cdr_positions[isotype]['CDR3'][1] - cdr_positions[isotype]['CDR3'][0],
    #                                      max(n) + 2, facecolor= cdr_color, lw = 0))
    #n, bins, p = pylab.hist(shm_pos, color = isotype_colors[isotype], bins = 50)
    plt.xlabel("Relative position of " + isotype + "V SHM in read", fontsize = 16)
    plt.ylabel("# SHMs", fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.title("SHMs in " + isotype + "V", fontsize = 18)
    # plot for SHM number
    plt.subplot(212)
    plt.hist(num_shms, color = isotype_colors[isotype], bins = 50, alpha = .75)
    plt.xlabel("# SHMs in " + isotype + "V", fontsize = 16)
    plt.ylabel("# sequences", fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    output_fname = output_prefix + "_" + isotype + "V"
    utils.output_figure(output_fname, "Distribution of # SHMs in " + isotype + "V segments", log)

def visualize_v_mutations_stats(shms_df, output_fname, log):
    all_shms_pos = {'IGH': [], 'IGK': [], 'IGL': []}
    num_all_shms = {'IGH': [], 'IGK': [], 'IGL': []}
    for it in shms_df:
        if not it.is_variable():
            continue
        read_shms = shms_df[it]
        for shm in read_shms:
            all_shms_pos[it.chain_type].append(float(shm.read_pos) / float(it.read_len))
        num_all_shms[it.chain_type].append(len(shms_df[it]))
    for iso in num_all_shms:
        if len(num_all_shms[iso]) < 10:
            log.info("# SHMs for " + iso + " is too small (" + str(len(num_all_shms[iso])) + "). Plot drawing was skipped")
            continue
        output_shm_stats_for_isotype(num_all_shms[iso], all_shms_pos[iso], iso, output_fname, log)

def get_aa_list():
    return ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def aa_is_valid(aa):
    return aa != '*' and aa != '-' and aa != 'X'

def get_aa_ticks_colors(aa_list):
    colors = []
    for aa in aa_list:
        if aa in utils.hydrophobic:
            colors.append('red')
        elif aa in utils.neutral:
            colors.append('green')
        elif aa in utils.hydrophilic:
            colors.append('blue')
    return colors

def visualize_aa_substitution_matrix(shms_df, output_fname, log):
    dict_aa = dict()
    num_shms = 0
    for it in shms_df:
        read_shms = shms_df[it]
        prev_pos = -1
        for shm in read_shms:
            if prev_pos / 3 != shm.read_pos / 3:
                if aa_is_valid(shm.gene_aa) and aa_is_valid(shm.read_aa):
                    aa_pair = shm.gene_aa + shm.read_aa
                    if not aa_pair in dict_aa:
                        dict_aa[aa_pair] = 0
                    dict_aa[aa_pair] += 1
                    num_shms += 1
            prev_pos = shm.read_pos
    aa_list = get_aa_list()
    aa_freq = []
    for i in range(0, len(aa_list)):
        aa_freq.append([0] * len(aa_list))
    for aa_pair in dict_aa:
        aa_freq[aa_list.index(aa_pair[1])][aa_list.index(aa_pair[0])] = float(dict_aa[aa_pair]) / float(num_shms)
    fig, ax = plt.subplots()
    sns.heatmap(aa_freq, cmap = plt.cm.jet, xticklabels = aa_list, yticklabels = aa_list, square = True, ax = ax)
    ax.tick_params(labelsize = 14)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12, rotation='horizontal')
    #for tick, color in zip(ax.get_xticklabels(), get_aa_ticks_colors(aa_list)):
    #    tick.set_color(color)
    #for tick, color in zip(ax.get_yticklabels(), get_aa_ticks_colors(aa_list)):
    #    tick.set_color(color)
    plt.xlabel("To", fontsize = 14)
    plt.ylabel("From", fontsize = 14, rotation='horizontal')
    utils.output_figure(output_fname, "Amino acid substitution heatmap", log)
    return aa_freq

def nucl_is_valid(nucl):
    return nucl != 'N'

def visualize_nucl_substitution_matrix(shms_df, output_fname, log):
    nucl_list = ['A', 'C', 'G', 'T']
    nucl_matrix = []
    for n in nucl_list:
        nucl_matrix.append([0] * len(nucl_list))
    num_shms = 0
    for it in shms_df:
        read_shms = shms_df[it]
        for shm in read_shms:
            if not shm.is_substitution():
                continue
            if nucl_is_valid(shm.read_nucl) and nucl_is_valid(shm.gene_nucl):
                nucl_matrix[nucl_list.index(shm.read_nucl)][nucl_list.index(shm.gene_nucl)] += 1
                num_shms += 1
    for i in range(0, len(nucl_matrix)):
        for j in range(0, len(nucl_matrix[i])):
            nucl_matrix[i][j] = float(nucl_matrix[i][j]) / float(num_shms)
    fig, ax = plt.subplots()
    sns.heatmap(nucl_matrix, cmap = plt.cm.Blues, xticklabels = nucl_list, yticklabels = nucl_list, square = True, ax = ax)
    ax.tick_params(labelsize = 14)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12, rotation='horizontal')
    plt.xlabel("To", fontsize = 14)
    plt.ylabel("From", fontsize = 14, rotation='horizontal')
    utils.output_figure(output_fname, "Nucleotide substitution heatmap", log)

def output_synonymous_shms(synonymous_pos, output_fname, log):
    if len(synonymous_pos) < 100:
        return
    plt.hist(synonymous_pos, color = 'r', bins = 100)
    plt.xlabel("Relative position of V SHM in read", fontsize = 14)
    plt.ylabel("#SHMs", fontsize = 14)
    plt.xlim(0, .75)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    utils.output_figure(output_fname, "Distribution of synonymous SHM positions in V segment", log)

def visualize_special_shm_positions(shm_df, syn_output_fname, special_output_fname, log):
    synonymous_pos = []
    stop_codon_pos = []
    deletion_pos = []
    insertion_pos = []
    for it in shm_df:
        read_shms = shm_df[it]
        for shm in read_shms:
            if not it.is_variable():
                continue
            relative_pos = float(shm.read_pos) / float(it.read_len)
            if shm.synonymous:
                synonymous_pos.append(relative_pos)
            elif shm.to_stop_codon:
                stop_codon_pos.append(relative_pos)
            elif shm.is_deletion():
                deletion_pos.append(relative_pos)
            elif shm.is_insertion():
                insertion_pos.append(relative_pos)
    output_synonymous_shms(synonymous_pos, syn_output_fname, log)
    pos = []
    labels = []
    colors = []
    plt.figure(figsize=(12, 9))
        #sns.distplot(synonymous_pos, hist = False, label = "Synonymous SHMs", color = 'r')
    #if len(stop_codon_pos) > 100:
    #    pos.append(stop_codon_pos)
    #    labels.append('Stop codon')
    #    colors.append('g')
    #    #sns.distplot(stop_codon_pos, hist = False, label = "Stop codon SHMs", color = 'g')
    if len(deletion_pos) > 10:
        pos.append(deletion_pos)
        labels.append('Deletions')
        colors.append('b')
        #sns.distplot(deletion_pos, hist = False, label = "Deletion SHMs", color = 'b')
    if len(insertion_pos) > 10:
        pos.append(insertion_pos)
        labels.append('Insertions')
        colors.append('g')
    if len(pos) == 0:
        log.info("Output contains very low number of special SHMs. Plot drawing will be skipped")
        return
    #sns.distplot(insertion_pos, hist = False, label = "Insertion SHMs", color = 'orange')
    plt.hist(pos, color = colors, label= labels, bins = 100 / len(pos))
    plt.xlim(0, .75)
    plt.legend(loc = 'upper center', ncol = len(pos), fontsize = 12, bbox_to_anchor=(0.5, -0.07))
    plt.xlabel("Relative position of V SHM in read", fontsize = 14)
    plt.ylabel("# SHMs", fontsize = 14)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    utils.output_figure(special_output_fname, "Distribution of indel V SHM positions in read", log)

def visualize_indel_shm_lengths(shm_df, output_fname, log):
    prev_read_pos = -1
    prev_gene_pos = -1
    insertion_length = []
    deletions_lengths = []
    in_len = 0
    del_len = 0
    for it in shm_df:
        read_shms = shm_df[it]
        for shm in read_shms:
            if shm.is_deletion():
                if shm.gene_pos - prev_gene_pos == 1:
                    del_len += 1
                else:
                    if del_len > 0:
                        deletions_lengths.append(del_len)
                    del_len = 1
                prev_gene_pos = shm.gene_pos
            if shm.is_insertion():
                if shm.read_pos - prev_read_pos == 1:
                    in_len += 1
                else:
                    if in_len > 0:
                        insertion_length.append(in_len)
                    in_len = 1
                prev_read_pos = shm.read_pos
    if in_len != 0:
        insertion_length.append(in_len)
    if del_len != 0:
        deletions_lengths.append(del_len)
    dt = []
    labels = []
    max_x_value = 0
    if len(deletions_lengths) > 10:
        dt.append(deletions_lengths)
        labels.append("Deletions")
    if len(insertion_length) > 10:
        dt.append(insertion_length)
        labels.append("Insertions")
    if len(dt) == 0:
        log.info("Output contains very low number of indel SHMs. Plot drawing was skipped")
        return
    plt.hist(dt, label = labels, bins = 50)
    plt.legend(loc = 'upper center', ncol = len(dt), fontsize = 14)
    plt.xlabel("Insertion / deletion SHM length", fontsize = 16)
    plt.ylabel("# insertion / deletion SHMs", fontsize = 16)
    xlim_right = 0
    if len(deletions_lengths) != 0:
        xlim_right = max(deletions_lengths)
    if len(insertion_length) != 0:
        xlim_right = max(xlim_right, max(insertion_length)) 
    plt.xlim(.5, xlim_right + .5)
    plt.xticks(range(0, xlim_right + 1), fontsize = 14)
    plt.yticks(fontsize = 14)
    utils.output_figure(output_fname, "Distribution of insertion/deletion SHM lengths", log)

def output_aa_freq(aa_freq, output_fname, log):
    aa_list = get_aa_list()
    fhandler = open(output_fname, "w")
    fhandler.write("from/to\t" + "\t".join(aa_list) + "\n")
    for i in range(0, len(aa_list)):
        fhandler.write(aa_list[i] + "\t" + "\t".join([str(ff) for ff in aa_freq[i]]) + "\n")
    log.info("Amino acid substitution matrix was written to " + output_fname)

def main(shm_df_fname, plot_dir, output_dir, log):
    log.info("== Output SHMs statistics")
    shm_df = SHMs(shm_df_fname)
    log.info(str(len(shm_df)) + " records were extracted from " + shm_df_fname)
    if len(shm_df) == 0:
        log.info("SHM data-frame contains 0 records. SHM visualization will be skipped")
        return
    visualize_v_mutations_stats(shm_df, os.path.join(plot_dir, "mutations_distribution"), log)
    aa_freq = visualize_aa_substitution_matrix(shm_df, os.path.join(plot_dir, "aa_substitutions"), log)
    output_aa_freq(aa_freq, os.path.join(output_dir, "aa_substitution_matrix.txt"), log)
    # todo: output aa freq
    visualize_nucl_substitution_matrix(shm_df, os.path.join(plot_dir, "nucl_substitutions"), log)
    visualize_special_shm_positions(shm_df, os.path.join(plot_dir, "synonymous_shms_positions"),
                                    os.path.join(plot_dir, "special_shms_positions"), log)
    visualize_indel_shm_lengths(shm_df, os.path.join(plot_dir, "indel_shms_length"), log)

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print "Invalid input parameters"
        print "python visualize_shm_stats.py shm_df.txt plot_dir output_dir logger"
        sys.exit(1)
    log = utils.get_logger_by_arg(sys.argv[4])
    main(sys.argv[1], sys.argv[2], sys.argv[3], log)
