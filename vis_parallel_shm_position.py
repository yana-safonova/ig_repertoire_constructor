import os
import sys

import matplotlib as mplt
mplt.use('Agg')

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class SHMRecord:
    def __init__(self, splits):
        self.edges = splits[8]
        self.codons = splits[9]

    def __eq__(self, other):
        return self.edges == other.edges and self.codons == other.codons

    def __hash__(self):
        return hash(self.edges) * hash(self.codons)

    def __str__(self):
        return self.edges + ": " + self.codons

    def __repr__(self):
        return str(self)

def ProcessFile(fname):
    lines = open(fname, "r").readlines()
    shms = dict()
    mult = dict()
    for i in range(2, len(lines)):
        splits = lines[i].strip().split()
        aa_pos = splits[1]
        aa = splits[2]
        record = SHMRecord(splits)
        pair = (aa_pos, aa)
        if pair not in shms:
            shms[pair] = set()
        shms[pair].add(record)
        if record not in mult:
            mult[record] = 0
        mult[record] += 1
    non_trivial_shms = 0
    non_trivial_aa = 0
    for s in shms:
        if len(shms[s]) > 1:
            print s, shms[s]
            non_trivial_aa += 1
            for r in shms[s]:
                non_trivial_shms += mult[r]
    print "# non-trivial SHMs: " + str(non_trivial_shms) + ", # supported AA: " + str(non_trivial_aa)

class SHM:
    def __init__(self, line):
        splits = line.strip().split()
        self.pos = int(splits[0])
        self.mult = int(splits[4])
        #self.trusted_mult = int(splits[2])
        self.hotspot = int(splits[5])
        self.stop_codon = int(splits[6])
        self.synonymous = int(splits[7])

    def Synonymous(self):
        return self.synonymous == 1

    def HotSpot(self):
        return self.hotspot == 1

def ProcessSHMFile(shm_file):
    lines = open(shm_file, "r").readlines()
    vgene_name = lines[0][len('#V_gene_name:'):]
    shms = []
    for i in range(2, len(lines)):
        shms.append(SHM(lines[i]))
    print str(len(lines) - 2) + " SHMs were extracted from " + shm_file
    return vgene_name, shms

def UpdateVgeneDict(v_gene, v_dict):
    v_gene = v_gene.split('*')[0]
    if v_gene not in v_dict:
        v_dict[v_gene] = 0
    v_dict[v_gene] += 1
    return v_gene, v_dict

def GetOutputFname(v_gene, v_dict, output_dir):
    return os.path.join(output_dir, v_gene + "_" + str(v_dict[v_gene]))

def OutputSHMList(shm_list, color, label):
    x = [s.pos for s in shm_list]
    y = [s.mult for s in shm_list]
    plt.plot(x, y, color=color, label=label, marker='.', linestyle='None', alpha = .75)

def OutputPlotInPdf(output_fname):
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()

def OutputSHMs(shms, output_fname):
    syn_shms = []
    hs_shms = []
    sc_shms = []
    other_shms = []
    max_mult = 0
    for s in shms:
        if s.Synonymous():
            syn_shms.append(s)
        elif s.HotSpot():
            hs_shms.append(s)
        elif s.stop_codon == 1:
            sc_shms.append(s)
        else:
            other_shms.append(s)
        max_mult = max(max_mult, s.mult)
    OutputSHMList(sc_shms, 'green', 'stop codons')
    OutputSHMList(syn_shms, "orange", 'synonymous')
    OutputSHMList(hs_shms, "blue", 'hot spots')
    OutputSHMList(other_shms, 'grey', 'other')
    plt.ylim(0, max_mult + 1)
    plt.xlabel('SHM position')
    plt.ylabel('SHM frequency')
    plt.legend(loc='upper center')
    OutputPlotInPdf(output_fname)

def OutputSHMPosPlot(shms, output_fname):
    shm_dict = dict()
    for s in shms:
        if s.pos not in shm_dict:
            shm_dict[s.pos] = 0
        shm_dict[s.pos] += s.mult
    shm_dict_sorted = sorted(shm_dict)
    x = shm_dict_sorted #[s for s in shm_dict]
    y = [shm_dict[s] for s in x] #[s for shm_dict[s] in shm_dict]
    plt.plot(x, y, color = 'blue', marker = '.')
    plt.plot(x, y, color = 'red', marker = '.', linestyle = 'None')
    OutputPlotInPdf(output_fname)

def FileTobeSkipped(shm_file):
    return shm_file[len(shm_file) - len('.DS_Store'):] == '.DS_Store' or \
           shm_file[len(shm_file) - len('.pdf'):] == '.pdf'

def main(shm_dir):
    shm_files = os.listdir(shm_dir)
    v_dict = dict()
    for f in shm_files:
        if FileTobeSkipped(f):
            continue
        print "== Processing " + f + "..."
        v_gene, shms = ProcessSHMFile(os.path.join(shm_dir, f))
        v_gene, v_dict = UpdateVgeneDict(v_gene, v_dict)
        output_base = GetOutputFname(v_gene, v_dict, shm_dir)
        OutputSHMs(shms, output_base + ".pdf")
        ProcessFile(os.path.join(shm_dir, f))
        print ""
        #OutputSHMPosPlot(shms, output_base + "_unique_shms.pdf")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print "Invalid input arguments"
        sys.exit(1)
    main(sys.argv[1])