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

# ----------------------- Structs for parsed output ---------------------
class RawIgblastOutput:
    blocks = list()
    read_index = dict()

    def GetBlockByName(self, query_name):
        if not query_name in self.read_index:
            print("ERROR: IgBlast output does not contain query " + query_name)
        return self.blocks[self.read_index[query_name]]

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

    from_ind = 1
    to_ind = 2
    len_ind = 3
    match_ind = 4
    mismatch_ind = 5
    gaps_ind = 6
    identity_ind = 7

    # general
    query_str = "# Query"
    db_str = "# Database"
    domain_str = "# Domain"
    vdj_rearrange_str = "# V-(D)-J rearrangement"
    vdj_juction_str = "# V-(D)-J junction"
    alignment_str = "# Alignment summary"
    hit_table_str = "# Hit table"

class VDJ_rearrangement:
    chain_type = ""
    top_v_genes = list()
    top_d_genes = list()
    top_j_genes = list()
    direct_strand = True

class IgRegionStats:
    from_ind = sys.maxint
    to_ind = -1
    length = 0
    matches = 0
    mismatches = 0
    gaps = 0
    percent_identity = 0
    valid = False

class FRs_CDRs:
    fr1 = IgRegionStats()
    cdr1 = IgRegionStats()
    fr2 = IgRegionStats()
    cdr2 = IgRegionStats()
    fr3 = IgRegionStats()
    cdr3 = IgRegionStats()

class BlockStats:
    query_name = ""
    db = list()
    domain = ""
    vdj_rearrangement = VDJ_rearrangement()
    vdj_junction = ""
    alignment_summary = FRs_CDRs()
    hit_table = list()

class BlockAuxStats:
    vdj_rearrange_ind = -1
    vdj_junction_ind = -1
    align_sum_start_ind = -1
    align_sum_end_ind = -1
    hit_table_start_ind = -1
    hit_table_end_ind = -1

# --------------------- Print Functions --------------------------

def PrintIgRegionStats(ig_region_stats):
    if not ig_region_stats.valid:
        print("Region is not found")
        return 
    string = "From: " + str(ig_region_stats.from_ind) + " " + "To: " + str(ig_region_stats.to_ind) + " " + "Length: " + str(ig_region_stats.length) + " " + "Matches: " + str(ig_region_stats.matches) + " " + "Mismatches: " + str(ig_region_stats.mismatches) +  " " + "Gaps: " + str(ig_region_stats.gaps) + " " + "% identity: " + str(ig_region_stats.percent_identity)
    print(string)

def PrintFRsCDRs(frs_cdrs):
    print("FR1:")
    PrintIgRegionStats(frs_cdrs.fr1)
    print("CDR1:")
    PrintIgRegionStats(frs_cdrs.cdr1)
    print("FR2:")
    PrintIgRegionStats(frs_cdrs.fr2)
    print("CDR2:")
    PrintIgRegionStats(frs_cdrs.cdr2)
    print("FR3:")
    PrintIgRegionStats(frs_cdrs.fr3)
    print("CDR3:")
    PrintIgRegionStats(frs_cdrs.cdr3)

def PrintVDJRearrangement(vdj_rearrangement):
    print("Chain type:" + vdj_rearrangement.chain_type)
    print("Top V genes:")
    print(vdj_rearrangement.top_v_genes)
    print("Top D genes:")
    print(vdj_rearrangement.top_d_genes)
    print("Top J genes:")
    print(vdj_rearrangement.top_j_genes)
    if vdj_rearrangement.direct_strand:
        print("Strand is direct")
    else:        
        print("Strand is complementary")

def PrintBlockStats(block_stats):
    print("Query name: " + block_stats.query_name)
    print("Database: ")
    print(block_stats.db)
    print("Domain: " + block_stats.domain)
    print("VDJ rearrangement: ")
    PrintVDJRearrangement(block_stats.vdj_rearrangement)
    print("VDJ junction: " + block_stats.vdj_junction)
    print("Alignment summary:")
    PrintFRsCDRs(block_stats.alignment_summary)
    #print("Hit table:")
    #print(block_stats.hit_table)

def PrintBlockAuxStats(block_aux_stats):
    print("vdj_rearrange_ind: " + str(block_aux_stats.vdj_rearrange_ind))
    print("vdj_junction_ind: " + str(block_aux_stats.vdj_junction_ind))
    print("alignment summary: " + str(block_aux_stats.align_sum_start_ind) + " " + str(block_aux_stats.align_sum_end_ind))
    print("hit table: " + str(block_aux_stats.hit_table_start_ind) + " " + str(block_aux_stats.hit_table_end_ind))

# -------------------------- Aux Functions --------------------------------

def LineStartMatchWithPrefix(line, prefix):
    return line[ : len(prefix)] == prefix

def CleanFRs_CDRs(frs_cdrs):
    frs_cdrs.fr1 = IgRegionStats()
    frs_cdrs.fr2 = IgRegionStats()
    frs_cdrs.fr3 = IgRegionStats()
    frs_cdrs.cdr1 = IgRegionStats()
    frs_cdrs.cdr2 = IgRegionStats()
    frs_cdrs.cdr3 = IgRegionStats()

# -------------------------- Line Process Functions -----------------------

def LineIsQuery(line):
    return LineStartMatchWithPrefix(line, BlockConfig.query_str) #line[ : len(BlockConfig.query_str)] == BlockConfig.query_str

def LineIsDatabase(line):
    return LineStartMatchWithPrefix(line, BlockConfig.db_str) #line[ : len(BlockConfig.db_str)] == BlockConfig.db_str

def LineIsDomain(line):
    return LineStartMatchWithPrefix(line, BlockConfig.domain_str) # line[ : len(BlockConfig.domain_str)] == BlockConfig.domain_str

def LineIsVDJRearrangeHeader(line):
    return LineStartMatchWithPrefix(line, BlockConfig.vdj_rearrange_str) #line[ : len(BlockConfig.vdj_rearrange_str)] == BlockConfig.vdj_rearrange_str

def LineIsVDJJuncHeader(line):
    return LineStartMatchWithPrefix(line, BlockConfig.vdj_juction_str) #line[ : len(BlockConfig.vdj_juction_str)] == BlockConfig.vdj_juction_str

def LineIsAlignSummaryHeader(line):
    return LineStartMatchWithPrefix(line, BlockConfig.alignment_str) #line[ : len(BlockConfig.alignment_str)] == BlockConfig.alignment_str

def LineIsHitTableHeader(line):
    return LineStartMatchWithPrefix(line, BlockConfig.hit_table_str) #line[ : len(BlockConfig.hit_table_str)] == BlockConfig.hit_table_str

# -------------

def ProcessQueryStr(line, block_stats):
    splits = line.split(' ')
    block_stats.query_name = splits[BlockConfig.query_name_ind].strip()

def ProcessDBStr(line, block_stats):
    splits = line.split(' ')
    for i in range(BlockConfig.db_start_ind, len(splits)):
        block_stats.db.append(splits[i].strip())

def ProcessDomainStr(line, block_stats):
    splits = line.split(' ')
    block_stats.domain = splits[BlockConfig.domain_ind].strip()

def VerifyBlockAuxStats(block_aux_stats):
    if block_aux_stats.align_sum_start_ind == -1:
        block_aux_stats.align_sum_end_ind = -1
    if block_aux_stats.align_sum_end_ind == -1:
        block_aux_stats.align_sum_start_ind = -1

    if block_aux_stats.hit_table_start_ind == -1:
        block_aux_stats.hit_table_end_ind = -1
    if block_aux_stats.hit_table_end_ind == -1:
        block_aux_stats.hit_table_start_ind = -1

# ------------------------- process VDJ rearrangement -------------------

def ProcessVDJRearrangemet(block, block_stats, block_aux_stats):
    if block_aux_stats.vdj_rearrange_ind == -1:
        return 
    vdj_rearrange_str = block[block_aux_stats.vdj_rearrange_ind].strip()

    block_stats.vdj_rearrangement = VDJ_rearrangement() 
    splits = vdj_rearrange_str.split('\t')

    block_stats.vdj_rearrangement.chain_type = splits[len(splits) - BlockConfig.chain_type_ind]

    if block_stats.vdj_rearrangement.chain_type == BlockConfig.heavy_chain_str:
        block_stats.vdj_rearrangement.top_v_genes = sorted(splits[BlockConfig.top_v_genes_hc_ind].strip().split(','))
        block_stats.vdj_rearrangement.top_d_genes = sorted(splits[BlockConfig.top_d_genes_hc_ind].strip().split(',')) 
        block_stats.vdj_rearrangement.top_j_genes = sorted(splits[BlockConfig.top_j_genes_hc_ind].strip().split(',')) 
    else:
        block_stats.vdj_rearrangement.top_v_genes = sorted(splits[BlockConfig.top_v_genes_lc_ind].strip().split(','))
        block_stats.vdj_rearrangement.top_j_genes = sorted(splits[BlockConfig.top_j_genes_lc_ind].strip().split(',')) 

    if splits[len(splits) - BlockConfig.strand_ind] == '-':
        block_stats.vdj_rearrangement.direct_strand = False

# ------------------------- process VDJ junction -------------------

def ProcessVDJJunction(block, block_stats, block_aux_stats):
    if block_aux_stats.vdj_junction_ind == -1:
        return 
    block_stats.vdj_junction = block[block_aux_stats.vdj_junction_ind].strip()

# ------------------------- process alignment summary ---------------

def ParseIgRegionStats(line):
    splits = line.split('\t')
    region = IgRegionStats()

    # from
    if splits[BlockConfig.from_ind] != "N/A":
        region.from_ind = int(splits[BlockConfig.from_ind])
    else:
        return region
    # to
    if splits[BlockConfig.to_ind] != "N/A":
        region.to_ind = int(splits[BlockConfig.to_ind])
    else:
        return region
    # length
    if splits[BlockConfig.len_ind] != "N/A":
        region.length = int(splits[BlockConfig.len_ind])
    else:
        return region
    # matches
    if splits[BlockConfig.match_ind] != "N/A":
        region.matches = int(splits[BlockConfig.match_ind])
    else:
        return region
    # mismatches
    if splits[BlockConfig.mismatch_ind] != "N/A":
        region.mismatches = int(splits[BlockConfig.mismatch_ind])
    else:
        return region
    # gaps
    if splits[BlockConfig.gaps_ind] != "N/A":
        region.gaps = int(splits[BlockConfig.gaps_ind])
    else:
        return region
    # % identity
    if splits[BlockConfig.identity_ind] != "N/A":
        region.percent_identity = float(splits[BlockConfig.identity_ind])
    else:
        return region

    region.valid = True
    return region

def ProcessAlignmentSummary(block, block_stats, block_aux_stats):
    if block_aux_stats.align_sum_start_ind == -1 or block_aux_stats.align_sum_end_ind == -1:
        return 

    frs_cdrs = FRs_CDRs()
    CleanFRs_CDRs(frs_cdrs)
    for i in range(block_aux_stats.align_sum_start_ind, block_aux_stats.align_sum_end_ind):
        line = block[i].strip()
        if LineStartMatchWithPrefix(line, BlockConfig.fr1_str):
            frs_cdrs.fr1 = ParseIgRegionStats(line)
        elif LineStartMatchWithPrefix(line, BlockConfig.fr2_str):
            frs_cdrs.fr2 = ParseIgRegionStats(line)
        elif LineStartMatchWithPrefix(line, BlockConfig.fr3_str):
            frs_cdrs.fr3 = ParseIgRegionStats(line)

        elif LineStartMatchWithPrefix(line, BlockConfig.cdr1_str):
            frs_cdrs.cdr1 = ParseIgRegionStats(line)
        elif LineStartMatchWithPrefix(line, BlockConfig.cdr2_str):
            frs_cdrs.cdr2 = ParseIgRegionStats(line)
        elif LineStartMatchWithPrefix(line, BlockConfig.cdr3_str):
            frs_cdrs.cdr3 = ParseIgRegionStats(line)
    block_stats.alignment_summary = frs_cdrs

# ------------------------- process hit table ----------------------

def ProcessHitTable(block, block_stats, block_aux_stats):
    if block_aux_stats.hit_table_start_ind == -1 or block_aux_stats.hit_table_end_ind == -1:
        return 
    for i in range(block_aux_stats.hit_table_start_ind, block_aux_stats.hit_table_end_ind + 1):
        block_stats.hit_table.append(block[i].strip())

# ------------------------- Read IgBlast --------------------------------

def ProcessBlock(block):
    if len(block) == 0:
        return BlockStats()

    block_stats = BlockStats()
    block_stats.db = list()
    block_aux_stats = BlockAuxStats()

    index = 0
    for l in block:
        if LineIsQuery(l):
            ProcessQueryStr(l, block_stats)
        elif LineIsDatabase(l):
            ProcessDBStr(l, block_stats)
        elif LineIsDomain(l):
            ProcessDomainStr(l, block_stats)
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
    VerifyBlockAuxStats(block_aux_stats)

    ProcessVDJRearrangemet(block, block_stats, block_aux_stats)
    ProcessVDJJunction(block, block_stats, block_aux_stats)
    ProcessAlignmentSummary(block, block_stats, block_aux_stats)
    ProcessHitTable(block, block_stats, block_aux_stats)

    #PrintBlockStats(block_stats)
    #PrintBlockAuxStats(block_aux_stats)
    #print("-----------")
    return block_stats

def CreateRawIgblastOutput(lines):
    igblast_output = RawIgblastOutput()
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

def ParseIgBlastOutput(igblast_fname, log):
    if not os.path.exists(igblast_fname):
        log.info("ERROR: IgBlast output file " + igblast_fname + " was not found")
    output = open(igblast_fname, "r")
    lines = output.readlines()
    igblast_output = CreateRawIgblastOutput(lines)
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

