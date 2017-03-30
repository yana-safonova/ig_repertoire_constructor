#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import getopt
import os
import logging
import shutil

from Bio.Seq import Seq, translate
from os.path import isfile, isdir, join
from os import listdir, curdir

from Bio import SeqIO
import subprocess

#########
class BlockConfig:
    query_name_ind = 2
    db_start_ind = 2
    domain_ind = 4

    # vdj rearrangement
    top_v_genes_hc_ind = 0
    top_d_genes_hc_ind = 1
    top_j_genes_hc_ind = 2
    top_v_genes_lc_ind = 0
    top_j_genes_lc_ind = 1
    # index from end of string
    chain_type_ind = 5
    strand_ind = 1

    heavy_chain_str = "VH"
    light_chain_str = "VL"

    # alignment summary
    fr1_str = "FR1"
    fr2_str = "FR2"
    fr3_str = "FR3"
    cdr1_str = "CDR1"
    cdr2_str = "CDR2"
    cdr3_str = "CDR3"

    from_ind = 7
    to_ind = 6
    len_ind = 5
    match_ind = 4
    mismatch_ind = 3
    gaps_ind = 2
    identity_ind = 1

    # general
    query_str = "# Query"
    db_str = "# Database"
    domain_str = "# Domain"
    vdj_rearrange_str = "# V-(D)-J rearrangement"
    vdj_juction_str = "# V-(D)-J junction"
    alignment_str = "# Alignment summary"
    hit_table_str = "# Hit table"

    # hit_table
    hit_type = 0
    hit_query_id = 1
    hit_subject_id = 2
    hit_perc_identity = 3
    hit_align_length = 4
    hit_mismatches = 5
    hit_gap_opens = 6
    hit_gaps = 7
    hit_q_start = 8
    hit_q_end = 9
    hit_s_start = 10
    hit_s_end = 11
    hit_evalue = 12
    hit_bit_score = 13
    hit_subject_length = 14

######### VDJ_rearrangement #########
class VDJ_rearrangement:
    def __init__(self):
        self.chain_type = ""
        self.top_v_genes = list()
        self.top_d_genes = list()
        self.top_j_genes = list()
        self.direct_strand = True

    def __str__(self):
        string = "Chain type:" + self.chain_type + ". Top V genes: " + self.top_v_genes
        strign += ", top D genes: " + self.top_d_genes + ", top J genes: " + self.top_j_genes
        if self.direct_strand:
            return string + ". Direct strand"
        return string + ". Reverse strand"

    def InitializeFromBlockLines(self, block_lines, block_aux_stats):
        if block_aux_stats.vdj_rearrange_ind == -1:
            return
        vdj_rearrange_str = block_lines[block_aux_stats.vdj_rearrange_ind].strip()
        splits = vdj_rearrange_str.split()
        self.chain_type = splits[len(splits) - BlockConfig.chain_type_ind]
        if self.chain_type == BlockConfig.heavy_chain_str:
            self.top_v_genes = sorted(splits[BlockConfig.top_v_genes_hc_ind].strip().split(','))
            self.top_d_genes = sorted(splits[BlockConfig.top_d_genes_hc_ind].strip().split(','))
            self.top_j_genes = sorted(splits[BlockConfig.top_j_genes_hc_ind].strip().split(','))
        else:
            self.top_v_genes = sorted(splits[BlockConfig.top_v_genes_lc_ind].strip().split(','))
            self.top_j_genes = sorted(splits[BlockConfig.top_j_genes_lc_ind].strip().split(','))
        if splits[len(splits) - BlockConfig.strand_ind] == '-':
            self.direct_strand = False

######### IgRegionStats #########
class IgRegionStats:
    def __init__(self):
        self.from_ind = sys.maxint
        self.to_ind = -1
        self.length = 0
        self.matches = 0
        self.mismatches = 0
        self.gaps = 0
        self.percent_identity = 0
        self.valid = False

    def __str__(self):
        if not self.valid:
            return "Region was not found"
        string = "From: " + str(self.from_ind) + ", to: " + str(self.to_ind) + ", length: " + str(self.length)
        string += ", matches: " + str(self.matches) + ", mismatches: " + str(self.mismatches) +  ", gaps: " + str(self.gaps)
        string += ", % identity: " + str(self.percent_identity)
        return string

    def InitializeFromString(self, string):
        splits = string.split()
        # from
        if splits[len(splits) - BlockConfig.from_ind] != "N/A":
            self.from_ind = int(splits[len(splits) - BlockConfig.from_ind]) - 1
        else:
            return
        # to
        if splits[len(splits) - BlockConfig.to_ind] != "N/A":
            self.to_ind = int(splits[len(splits) - BlockConfig.to_ind]) - 1
        else:
            return
        # length
        if splits[len(splits) - BlockConfig.len_ind] != "N/A":
            self.length = int(splits[len(splits) - BlockConfig.len_ind])
        else:
            return
        # matches
        if splits[len(splits) - BlockConfig.match_ind] != "N/A":
            self.matches = int(splits[len(splits) - BlockConfig.match_ind])
        else:
            return
        # mismatches
        if splits[len(splits) - BlockConfig.mismatch_ind] != "N/A":
            self.mismatches = int(splits[len(splits) - BlockConfig.mismatch_ind])
        else:
            return
        # gaps
        if splits[len(splits) - BlockConfig.gaps_ind] != "N/A":
            self.gaps = int(splits[len(splits) - BlockConfig.gaps_ind])
        else:
            return
        # % identity
        if splits[len(splits) - BlockConfig.identity_ind] != "N/A":
            self.percent_identity = float(splits[len(splits) - BlockConfig.identity_ind])
        else:
            return
        self.valid = True

######### FRs_CDRs #########
class FRs_CDRs:
    def __init__(self):
        self.fr1 = IgRegionStats()
        self.cdr1 = IgRegionStats()
        self.fr2 = IgRegionStats()
        self.cdr2 = IgRegionStats()
        self.fr3 = IgRegionStats()
        self.cdr3 = IgRegionStats()

    def __str__(self):
        string = "FR1: " + str(self.fr1) + "\n"
        string += "CDR1: " + str(self.cdr1) + "\n"
        string += "FR2: " + str(self.fr2) + "\n"
        string += "CDR2: " + str(self.cdr2) + "\n"
        string += "FR3: " + str(self.fr3) + "\n"
        return string + "CDR3: " + str(self.cdr3)

    def InitializeFromBlockLines(self, block_lines, block_aux_stats):
        if block_aux_stats.align_sum_start_ind == -1 or block_aux_stats.align_sum_end_ind == -1:
            return
        for i in range(block_aux_stats.align_sum_start_ind, block_aux_stats.align_sum_end_ind):
            line = block_lines[i].strip()
            if LineStartMatchWithPrefix(line, BlockConfig.fr1_str):
                self.fr1.InitializeFromString(line)
            elif LineStartMatchWithPrefix(line, BlockConfig.fr2_str):
                self.fr2.InitializeFromString(line)
            elif LineStartMatchWithPrefix(line, BlockConfig.fr3_str):
                self.fr3.InitializeFromString(line)
            elif LineStartMatchWithPrefix(line, BlockConfig.cdr1_str):
                self.cdr1.InitializeFromString(line)
            elif LineStartMatchWithPrefix(line, BlockConfig.cdr2_str):
                self.cdr2.InitializeFromString(line)
            elif LineStartMatchWithPrefix(line, BlockConfig.cdr3_str):
                self.cdr3.InitializeFromString(line)

######### BlockAuxStats #########
class BlockAuxStats:
    def __init__(self):
        self.vdj_rearrange_ind = -1
        self.vdj_junction_ind = -1
        self.align_sum_start_ind = -1
        self.align_sum_end_ind = -1
        self.hit_table_start_ind = -1
        self.hit_table_end_ind = -1

    def __str__(self):
        return "VDJ_rearrange_ind: " + str(self.vdj_rearrange_ind) + ", VDJ_junction_ind: " + str(self.vdj_junction_ind) + ". Alignment summary: " + str(self.align_sum_start_ind) + " - " + str(self.align_sum_end_ind) + ". Hit table: " + str(self.hit_table_start_ind) + " - " + str(self.hit_table_end_ind)

    def Verify(self):
        if self.align_sum_start_ind == -1:
            self.align_sum_end_ind = -1
        if self.align_sum_end_ind == -1:
            self.align_sum_start_ind = -1
        if self.hit_table_start_ind == -1:
            self.hit_table_end_ind = -1
        if self.hit_table_end_ind == -1:
            self.hit_table_start_ind = -1

######### HitTableRow #########
class HitTableRow:
    def __init__(self, hit_string):
        splits = hit_string.split()
        self.type = splits[BlockConfig.hit_type]
        self.query_id = splits[BlockConfig.hit_query_id]
        self.subject_id = splits[BlockConfig.hit_subject_id]
        self.perc_identity = float(splits[BlockConfig.hit_perc_identity])
        self.align_length = int(splits[BlockConfig.hit_align_length])
        self.mismatches = int(splits[BlockConfig.hit_mismatches])
        self.gap_open = int(splits[BlockConfig.hit_gap_opens])
        self.gaps = int(splits[BlockConfig.hit_gaps])
        self.q_start = int(splits[BlockConfig.hit_q_start]) - 1
        self.q_end = int(splits[BlockConfig.hit_q_end]) - 1
        self.s_start = int(splits[BlockConfig.hit_s_start]) - 1
        self.s_end = int(splits[BlockConfig.hit_s_end]) - 1
        self.evalue = float(splits[BlockConfig.hit_evalue])
        self.bit_score = float(splits[BlockConfig.hit_bit_score])

    def __str__(self):
        return self.type + ": " + self.query_id + ", " + self.subject_id + ", " + str(self.perc_identity) + ", " + str(self.align_length) + ", " + str(self.mismatches) + ", " + str(self.gap_open) + ", " + str(self.gaps) + ", " + str(self.q_start) + ", " + str(self.q_end) + ", " + str(self.s_start) + ", " + str(self.s_end) + ", " + str(self.evalue) + ", " + str(self.bit_score)

class HitTable:
    def __init__(self):
        self.rows = list()
        self.index = -1

    def AddRow(self, hit_string):
        self.rows.append(HitTableRow(hit_string))

    def __str__(self):
        res = ""
        for row in self.rows:
            res += str(row) + "\n"
        return res

    def __len__(self):
        return len(self.rows)

    def __iter__(self):
        return self

    def next(self):
        if self.index + 1 < len(self):
            self.index += 1
            return self.rows[self.index]
        else:
            raise StopIteration

    def InitializeFromBlockLines(self, block_lines, block_aux_stats):
        if block_aux_stats.hit_table_start_ind == -1 or block_aux_stats.hit_table_end_ind == -1:
            return
        for i in range(block_aux_stats.hit_table_start_ind, block_aux_stats.hit_table_end_ind + 1):
            if block_lines[i][0] == "#":
                return
            self.AddRow(block_lines[i].strip())

######### IgblastBlock #########
class IgblastBlock:
    def __init__(self):
        self.query_name = ""
        self.db = list()
        self.domain = ""
        self.vdj_rearrangement = VDJ_rearrangement()
        self.vdj_junction = ""
        self.alignment_summary = FRs_CDRs()
        self.hit_table = HitTable()

    def __str__(self):
        string = "Query name: " + self.query_name + ", database: " + self.db + ", domain: " + self.domain + "\n"
        string += "VDJ rearrangement: " + str(self.vdj_rearrangement) + ". VDJ junction: " + str(self.vdj_junction) + "\n"
        string += "Alignment summary:" + str(self.alignment_summary) + "\n"
        string += "Hit table:" + str(self.hit_table)
        return string

    def ProcessQueryStr(self, line):
        splits = line.split(' ')
        self.query_name = splits[BlockConfig.query_name_ind].strip()

    def ProcessDBStr(self, line):
        splits = line.split(' ')
        for i in range(BlockConfig.db_start_ind, len(splits)):
            self.db.append(splits[i].strip())

    def ProcessDomainStr(self, line):
        splits = line.split(' ')
        self.domain = splits[BlockConfig.domain_ind].strip()

    def ProcessVDJRearrangemet(self, block_lines, block_aux_stats):
        self.vdj_rearrangement.InitializeFromBlockLines(block_lines, block_aux_stats)

    def ProcessVDJJunction(self, block_lines, block_aux_stats):
        if block_aux_stats.vdj_junction_ind == -1:
            return
        self.vdj_junction = block_lines[block_aux_stats.vdj_junction_ind].strip()

    def ProcessAlignmentSummary(self, block_lines, block_aux_stats):
        self.alignment_summary.InitializeFromBlockLines(block_lines, block_aux_stats)

    def ProcessHitTable(self, block_lines, block_aux_stats):
        self.hit_table.InitializeFromBlockLines(block_lines, block_aux_stats)

######### IgblastOutput #########
class IgblastOutput:
    def __init__(self):
        self.blocks = list()
        self.read_index = dict()

    def __getitem__(self, query_name):
        if not query_name in self.read_index:
            print("ERROR: IgBlast output does not contain query " + query_name)
            return
        return self.blocks[self.read_index[query_name]]

# -------------------------- Aux Functions --------------------------------
def LineStartMatchWithPrefix(line, prefix):
    return line[ : len(prefix)] == prefix

def LineIsQuery(line):
    return LineStartMatchWithPrefix(line, BlockConfig.query_str)

def LineIsDatabase(line):
    return LineStartMatchWithPrefix(line, BlockConfig.db_str)

def LineIsDomain(line):
    return LineStartMatchWithPrefix(line, BlockConfig.domain_str)

def LineIsVDJRearrangeHeader(line):
    return LineStartMatchWithPrefix(line, BlockConfig.vdj_rearrange_str)

def LineIsVDJJuncHeader(line):
    return LineStartMatchWithPrefix(line, BlockConfig.vdj_juction_str)

def LineIsAlignSummaryHeader(line):
    return LineStartMatchWithPrefix(line, BlockConfig.alignment_str)

def LineIsHitTableHeader(line):
    return LineStartMatchWithPrefix(line, BlockConfig.hit_table_str)

# ------------------------- Read IgBlast --------------------------------

def ProcessBlock(block):
    if len(block) == 0:
        return IgblastBlock()

    block_stats = IgblastBlock()
    block_aux_stats = BlockAuxStats()

    index = 0
    for l in block:
        if LineIsQuery(l):
            block_stats.ProcessQueryStr(l)
        elif LineIsDatabase(l):
            block_stats.ProcessDBStr(l)
        elif LineIsDomain(l):
            block_stats.ProcessDomainStr(l)
        elif LineIsVDJRearrangeHeader(l):
            block_aux_stats.vdj_rearrange_ind = index + 1
        elif LineIsVDJJuncHeader(l):
            block_aux_stats.vdj_junction_ind = index + 1
        elif LineIsAlignSummaryHeader(l):
            block_aux_stats.align_sum_start_ind = index + 1
        elif LineIsHitTableHeader(l):
            block_aux_stats.align_sum_end_ind = index - 2
            block_aux_stats.hit_table_start_ind = index + 3
        index += 1
    block_aux_stats.hit_table_end_ind = len(block) - 1
    block_aux_stats.Verify()

    block_stats.ProcessVDJRearrangemet(block, block_aux_stats)
    block_stats.ProcessVDJJunction(block, block_aux_stats)
    block_stats.ProcessAlignmentSummary(block, block_aux_stats)
    block_stats.ProcessHitTable(block, block_aux_stats)

    #PrintBlockStats(block_stats)
    #PrintBlockAuxStats(block_aux_stats)
    #print("-----------")
    return block_stats

def CreateIgblastOutput(lines):
    igblast_output = IgblastOutput()
    start_block_line = "# IGBLASTN"
    end_output_line = "# BLAST processed"
    block = list()
    for l in lines:
        if l[:len(start_block_line)] == start_block_line and len(block) != 0:
            block_stats = ProcessBlock(block)
            igblast_output.blocks.append(block_stats)
            block = list()
        else:
            block.append(l)
        if l[:len(end_output_line)] == end_output_line:
            break
    block_stats = ProcessBlock(block)
    igblast_output.blocks.append(block_stats)
    return igblast_output

def CreateReadIndexMap(igblast_output):
    index = 0
    for b in igblast_output.blocks:
        #print(b)
        #PrintBlockStats(b)
        igblast_output.read_index[b.query_name] = index
        index += 1
    #print(igblast_output.read_index)

def ParseIgBlastOutput(igblast_fname, log, openfnc=open):
    if not os.path.exists(igblast_fname):
        log.info("ERROR: IgBlast output file " + igblast_fname + " was not found")
    output = openfnc(igblast_fname, "r")
    lines = output.readlines()
    igblast_output = CreateIgblastOutput(lines)
    CreateReadIndexMap(igblast_output)
    log.info(str(len(igblast_output.blocks)) + " alignment block(s) were read from " + igblast_fname)
    log.info(str(len(igblast_output.read_index)) + " reads were processed")
    return igblast_output

# ---------------------- V(J) groups -----------------------------------

# temporary solution
def GeneNameWithoutVariations(name):
    name = name.replace('/', '-')
    splits = name.split('*')
    return splits[0]

def AddEdgeToTopGenesGraph(graph, gene1, gene2):
    if gene1 not in graph:
        graph[gene1] = {}
    if gene2 not in graph[gene1]:
        graph[gene1][gene2] = 0
    graph[gene1][gene2] += 1

def BuildTopGenesGraph(igblast_output, gene_type):
    graph = {}

    for block in igblast_output.blocks:
        top_genes = []
        if gene_type == 'V':
            top_genes = list(set([GeneNameWithoutVariations(name) for name in block.vdj_rearrangement.top_v_genes]))
        else:
            top_genes = list(set([GeneNameWithoutVariations(name) for name in block.vdj_rearrangement.top_j_genes]))
        if len(top_genes) == 1 and top_genes[0] not in graph:
            graph[top_genes[0]] = {}
        for i in xrange(len(top_genes)):
            for j in xrange(i + 1, len(top_genes)):
                AddEdgeToTopGenesGraph(graph, top_genes[i], top_genes[j])
                AddEdgeToTopGenesGraph(graph, top_genes[j], top_genes[i])

    return graph


def GroupVJGenes(igblast_output):
    align_groups = {}

    connections = {}
    connections['V'] = BuildTopGenesGraph(igblast_output, 'V')
    connections['J'] = BuildTopGenesGraph(igblast_output, 'J')

    weight_threshold = 10

    last_group_id = 0
    for gene_type, graph in connections.items():
        for gene, adj_list in graph.items():
            if gene in align_groups:
                continue
            last_group_id += 1
            align_groups[gene] = last_group_id
            queue = [g for g, w in adj_list.items() if w >= weight_threshold]
            i = 0
            while i != len(queue):
                if queue[i] in align_groups:
                    i += 1
                    continue
                align_groups[queue[i]] = last_group_id
                queue += [g for g, w in graph[queue[i]].items() if w >= weight_threshold]
                i += 1

    return align_groups

def VJGroupName(vgene, jgene):
    return vgene + "_" + jgene

def GetVGroupForRead(igblast_output, genes_align_groups, read_name, log):
    if not read_name in igblast_output.read_index:
        log.info("ERROR: Read " + read_name + " was not found in IgBlast output")
        sys.exit(1)
    index = igblast_output.read_index[read_name]
    block = igblast_output.blocks[index]

    if len(block.vdj_rearrangement.top_v_genes) == 0:
        log.info("V genes for read " + read_name + " were not determined")
        sys.exit(1)
    v_best_gene = block.vdj_rearrangement.top_v_genes[0]
    v_groups_name = "V" + str(genes_align_groups[GeneNameWithoutVariations(v_best_gene)])
    return v_groups_name

def GetVJGroupForRead(igblast_output, genes_align_groups, read_name, log):
    if not read_name in igblast_output.read_index:
        log.info("ERROR: Read " + read_name + " was not found in IgBlast output")
        sys.exit(1)
    index = igblast_output.read_index[read_name]
    block = igblast_output.blocks[index]

    if len(block.vdj_rearrangement.top_v_genes) == 0:
        log.info("V genes for read " + read_name + " were not determined")
        sys.exit(1)
    v_best_gene = block.vdj_rearrangement.top_v_genes[0]
    v_groups_name = "V" + str(genes_align_groups[GeneNameWithoutVariations(v_best_gene)])

    if len(block.vdj_rearrangement.top_j_genes) == 0:
        log.info("J genes for read " + read_name + " were not determined")
        sys.exit(1)
    j_best_gene = block.vdj_rearrangement.top_j_genes[0]
    j_groups_name = "J" + str(genes_align_groups[GeneNameWithoutVariations(j_best_gene)])
    return VJGroupName(v_groups_name, j_groups_name)

