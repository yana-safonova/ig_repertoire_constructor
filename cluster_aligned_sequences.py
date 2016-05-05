import os
import sys
import shutil

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import operator

import numpy as np
import scipy.cluster.hierarchy as hac
#import matplotlib.pyplot as plt

def ComputeOverlap(alignment, row1, row2):
    start_index = 0
    num_columns = len(alignment[0])
    end_index = num_columns - 1
    for i in range(0, num_columns):
        if alignment[row1][i] == '-' or alignment[row2][i] == '-':
            start_index = i
        else:
            break
    for i in range(0, num_columns):
        j = num_columns - i - 1
        if alignment[row1][j] == '-' or alignment[row2][j] == '-':
            end_index = j
        else:
            break
    return start_index, end_index


def ComputeDistance(alignment, row1, row2):
    distance = 0
    start, end = ComputeOverlap(alignment, row1, row2)
    for i in range(start, end + 1):
        if alignment[row1][i] != alignment[row2][i] and alignment[row1][i] != "N" and alignment[row2][i] != "N":
            distance += 1
    return distance

def ComputeDistanceMatrix(alignment):
    num_columns = len(alignment[0])
    num_rows = len(alignment)
    dist_matrix = [0] * (num_rows * (num_rows - 1) / 2)
    for i in range(0, (num_rows - 1)):
        for j in range(i + 1, num_rows):
            dist = ComputeDistance(alignment, i, j)
            shift = j - 1 - i
            dist_matrix[(2 * num_rows - 1 - i) * i / 2 + shift] = dist
    return dist_matrix

def ComputeSimpleDecompositions(link_matrix, threshold):
    clusters = [0] * (len(link_matrix) + 1)
    for i in range(0, len(clusters)):
        clusters[i] = i
    for i in range(0, len(link_matrix)):
        if link_matrix[i][2] > threshold:
            break
        ind1 = int(link_matrix[i][0])
        ind2 = int(link_matrix[i][1])
        new_cluster = len(clusters)
        clusters.append(new_cluster)
        clusters[ind1] = new_cluster
        clusters[ind2] = new_cluster
    clusters_dict = dict()
    for i in range(0, len(link_matrix) + 1):
        root = clusters[i]
        while clusters[root] != root:
            root = clusters[root]
        if root not in clusters_dict:
            clusters_dict[root] = set()
        clusters_dict[root].add(i)
    return clusters_dict

def ComputeSymbolDict(alignment, indices, pos):
    sym_dict = dict()
    for ind in indices:
        curr_symbol = alignment[ind][pos]
        if curr_symbol not in sym_dict:
            sym_dict[curr_symbol] = 0
        sym_dict[curr_symbol] += 1
    return sym_dict

def SymbolIsValid(sorted_dict, pos):
    return sorted_dict[pos][0] != '-' and sorted_dict[pos][0] != 'N'

def SymbolIsTrusted(sorted_dict, pos, threshold):
    return sorted_dict[pos][1] >= threshold

def ComputeConsensusForCluster(alignment, indices):
    consensus = ""
    trusted_threshold = min(len(indices), 5)
    for i in range(0, len(alignment[0])):
        sym_dict = ComputeSymbolDict(alignment, indices, i)
        sorted_dict = sorted(sym_dict.items(), key=operator.itemgetter(1), reverse=True)
        #print sorted_dict
        if SymbolIsValid(sorted_dict, 0):
            consensus += sorted_dict[0][0]
        else:
            for j in range(1, len(sorted_dict)):
                if SymbolIsValid(sorted_dict, j) and SymbolIsTrusted(sorted_dict, j, trusted_threshold):
                    consensus += sorted_dict[j][0]
                    break
    #print consensus
    return consensus

class Consensus:
    def _ComputeIsotypes(self, record_ids, indices):
        isotypes = set()
        for ind in indices:
            splits = record_ids[ind].split("|")
            splits = splits[1].split('_')
            if len(splits) == 1:
                splits = splits[0].split(':')
            isotypes.add(splits[1])
        return isotypes

    def _ComputeMBs(self, record_ids, indices):
        mbs = list()
        for ind in indices:
            splits = record_ids[ind].split("|")
            splits = splits[0].split('_')
            if len(splits) == 1:
                splits = splits[0].split(':')
            mbs.append(splits[1])
        return mbs

    def __init__(self, ind, seq, record_ids, indices):
        self.index = ind
        self.seq = str(seq)
        self.count = len(indices)
        self.isotypes = self._ComputeIsotypes(record_ids, indices)
        self.isotypes_str = ",".join(self.isotypes)
        self.mol_barcodes = self._ComputeMBs(record_ids, indices)

    def consensus_id(self):
        return "ID=" + str(self.index) + "|ISOTYPE=" + self.isotypes_str + "|COUNT=" + str(self.count)

    def __str__(self):
        return self.consensus_id() + "\n" + self.seq

def ComputeConsensusSequences(link_matrix, alignment):
    consensus_list = list()
    threshold = 10
    clusters = ComputeSimpleDecompositions(link_matrix, threshold)
    print "  " + str(len(clusters)) + " cluster(s) were constructed using threshold " + str(threshold)
    record_ids = list()
    for record in alignment:
        record_ids.append(record.id)
    index = 1
    for cluster_id in clusters:
        #print clusters[cluster_id]
        consensus_seq = ComputeConsensusForCluster(alignment, clusters[cluster_id])
        consensus = Consensus(index, consensus_seq, record_ids, clusters[cluster_id])
        consensus_list.append(consensus)
        #print consensus
        index += 1
    return consensus_list
    #print clusters

def PerformClusteringByAlignment(alignment):
    dist_matrix = ComputeDistanceMatrix(alignment)
    link_matrix = hac.single(dist_matrix)
    return ComputeConsensusSequences(link_matrix, alignment)

def WriteConsensusList(output_dir, prefix, consensus_list):
    splits = prefix.split("_")
    output_prefix = splits[0] + "_" + splits[1] + "_" + splits[2] + "_seqs_" + str(len(consensus_list))
    output_fasta = os.path.join(output_dir, output_prefix + ".fasta")
    output_map = os.path.join(output_dir, output_prefix + ".map")
    print "  Output of consensus sequences in " + output_fasta
    print "  Output of map for molecular barcodes in " + output_map
    fasta_handler = open(output_fasta, "w")
    map_handler = open(output_map, "w")
    for consensus in consensus_list:
        SeqIO.write(SeqRecord(Seq(consensus.seq), id = consensus.consensus_id(), description=""), fasta_handler, "fasta")
        for mb in consensus.mol_barcodes:
            map_handler.write(mb + "\t" + str(consensus.index) + "\n")
    fasta_handler.close()
    map_handler.close()

def main(argv):
    input_dir = argv[1]
    file_list = os.listdir(input_dir)
    prefix_set = set()
    for f in file_list:
        prefix_set.add(f.split(".")[0])
    print "Input dir " + input_dir + " contains " + str(len(prefix_set)) + " chain groups"

    output_dir = argv[2]
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    print "Output dir: " + output_dir
    os.mkdir(output_dir)
    for p in prefix_set:
        print "Processing files " + p + "..."
        fasta_file = os.path.join(input_dir, p + ".fasta")
        records = list(SeqIO.parse(open(fasta_file, "r"), "fasta"))
        consensus_list = list()
        if len(records) == 1:
            consensus_list.append(Consensus(1, records[0].seq, [records[0].id], {0}))
        else:
            align_file = os.path.join(input_dir, p + ".aln")
            alignment = AlignIO.read(align_file, "clustal")
            consensus_list = PerformClusteringByAlignment(alignment)
        WriteConsensusList(output_dir, p, consensus_list)

if __name__ == "__main__":
    main(sys.argv)
