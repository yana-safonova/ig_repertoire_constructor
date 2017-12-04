import os
import sys
import shutil
import matplotlib as mplt
#mplt.use('Agg')
import pandas as pd

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.patches as patches

class SHM:
    def __init__(self, line):
        splits = line.strip().split()
        self.pos = int(splits[0])
        self.aa_pos = int(splits[1])
        self.src_nucl = splits[2]
        self.dst_nucl = splits[3]
        self.src_aa = splits[4]
        self.dst_aa = splits[5]
        self.src_triplet = splits[6]
        self.dst_triplet = splits[7]
        self.mult = int(splits[8])
        self.region = splits[9]

    def synonymous(self):
        return self.src_aa == self.dst_aa

    def stop_codon(self):
        return self.dst_aa == '*'

    def __repr__(self):
        return self.dst_aa

    def __str__(self):
        return self.dst_aa

class Tree:
    def _init_cdrs(self, cdr1_line, cdr2_line, cdr3_line):
        self.cdr1_start = int(cdr1_line.split(',')[0][len('@CDR1:'):])
        self.cdr1_end = int(cdr1_line.split(',')[1])
        self.cdr2_start = int(cdr2_line.split(',')[0][len('@CDR2:'):])
        self.cdr2_end = int(cdr2_line.split(',')[1])
        self.cdr3_start = int(cdr3_line.split(',')[0][len('@CDR3:'):])
        self.cdr3_end = int(cdr3_line.split(',')[1])

    def _init_len(self, len_line):
        self.len = int(len_line[len('@VDJ_length:'):])

    def __init__(self, tree_fname):
        lines = open(tree_fname, "r").readlines()
        self._init_cdrs(lines[0], lines[1], lines[2])
        self._init_len(lines[3])
        self.shms = []
        for i in range(5, len(lines)):
            shm = SHM(lines[i])
            self.shms.append(shm)

    def max_multiplicity(self):
        max_mult = 0
        for s in self.shms:
            max_mult = max(max_mult, s.mult)
        return max_mult

def output_tree_shms(tree_fname, output_base):
    tree = Tree(tree_fname)
    sc_shms = [s for s in tree.shms if s.stop_codon()]
    syn_shms = [s for s in tree.shms if s.synonymous()]
    f, ax = plt.subplots()
#    plt.plot([s.pos for s in tree.shms], [s.mult for s in tree.shms], color = 'blue', marker = '.', linestyle = 'None', label = 'other')
#    plt.plot([s.pos for s in syn_shms], [s.mult for s in syn_shms], color = 'green', linestyle = 'None', marker = '.', label = 'synonymous')
#    plt.plot([s.pos for s in sc_shms], [s.mult for s in sc_shms], color = 'orange', marker = '.', linestyle = 'None', label = 'stop codon')

    x_dict = {'syn' : [], 'sc' : [], 'other' : []}
    y_dict = {'syn' : [], 'sc' : [], 'other' : []}
    for shm in tree.shms:
        if shm.synonymous():
            x_dict['syn'].append(shm.pos)
            y_dict['syn'].append(shm.mult)
        elif shm.stop_codon():
            x_dict['sc'].append(shm.pos)
            y_dict['sc'].append(shm.mult)
        else:
            x_dict['other'].append(shm.pos)
            y_dict['other'].append(shm.mult)

    plt.plot(x_dict['sc'], y_dict['sc'], color = 'orange', marker = '.', linestyle = 'None', label = 'stop codon')
    plt.plot(x_dict['syn'], y_dict['syn'], color = 'green', linestyle = 'None', marker = '.', label = 'synonymous')
    plt.plot(x_dict['other'], y_dict['other'], color = 'blue', marker = '.', linestyle = 'None', label = 'other')

    #aa_x_dict = dict()
    #aa_y_dict = dict()
    #for shm in tree.shms:
    #    aa = shm.dst_aa
    #    if aa not in aa_x_dict:
    #        aa_x_dict[aa] = []
    #        aa_y_dict[aa] = []
    #    aa_x_dict[aa].append(shm.pos)
    #    aa_y_dict[aa].append(shm.mult)

    #cmap = plt.cm.get_cmap('jet', 21)
    #i = 0
    #for aa in aa_x_dict:
    #    plt.plot(aa_x_dict[aa], aa_y_dict[aa], marker = '.', linestyle = 'None', label = aa, color = cmap(i))
    #    i += 1

    add_mult = int(float(tree.max_multiplicity()) / 100.0 * 10)
    ax.add_patch(patches.Rectangle((tree.cdr1_start, 0), tree.cdr1_end - tree.cdr1_start + 1, tree.max_multiplicity() + add_mult, alpha = 0.2, facecolor="red"))
    ax.add_patch(patches.Rectangle((tree.cdr2_start, 0), tree.cdr2_end - tree.cdr2_start + 1, tree.max_multiplicity() + add_mult, alpha = 0.2, facecolor="red"))
    ax.add_patch(patches.Rectangle((tree.cdr3_start, 0), tree.cdr3_end - tree.cdr3_start + 1, tree.max_multiplicity() + add_mult, alpha = 0.2, facecolor="red"))
    plt.xlabel('SHM position')
    plt.ylabel('SHM multiplicity')
    plt.xlim(0, tree.len + 5)
    plt.ylim(0, tree.max_multiplicity() + add_mult)
    plt.legend(loc='upper center', ncol = 3, fontsize='10', facecolor = "white", framealpha = 1, frameon=True)
    pp = PdfPages(output_base + "_shm_pos.pdf")
    pp.savefig()
    pp.close()
    plt.clf()
    print "Output of SHM plot for " + tree_fname
    
    df = pd.DataFrame({'mult': [s.mult for s in tree.shms], 'reg': [s.region for s in tree.shms], 'syn' : [s.synonymous() for s in tree.shms], 'aa' : [s.dst_aa for s in tree.shms]})
    sns.violinplot(x = 'reg', y = 'mult', data = df, order=["FR1", "CDR1", 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4'], linewidth = 0.5, palette=sns.color_palette("RdYlBu"))
    #sns.violinplot(x = 'aa', y = 'mult', data = df)
    plt.xlabel('Region')
    plt.ylabel('SHM multiplicity')
    pp = PdfPages(output_base + "_violin.pdf")
    pp.savefig()
    pp.close()
    plt.clf()

def GetTreeSize(tree_fname):
    return int(tree_fname.split('_')[4])

def main(input_dir, output_dir):
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)
    files = os.listdir(input_dir)
    for f in files:
        tree_fname = os.path.join(input_dir, f)
        if GetTreeSize(f) < 10:
            continue
        output_fname = os.path.join(output_dir, f.split('.')[0])
        output_tree_shms(tree_fname, output_fname)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
