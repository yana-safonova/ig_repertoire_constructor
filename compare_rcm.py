import os
import sys
import shutil

graph_utils = "src/graph_utils"
sys.path.append(graph_utils)

import graph_structs
from graph_structs import Graph
import graph_io

class Repertoire:
    def __ExtractFromRCM(self):
        fhandler = open(self.__rcm_fname)
        for line in fhandler:
            splits = line.strip().split()
            if len(splits) != 2:
                print "RCM file " + self.__rcm_fname + " is incorrect"
            read_name = splits[0]
            cluster_id = int(splits[1])
            if not cluster_id in self.__clusters_dict:
                self.__clusters_dict[cluster_id] = list()
            self.__clusters_dict[cluster_id].append(read_name)
            self.__read_cluster_map[read_name] = cluster_id
        fhandler.close()

    def __init__(self, rcm_fname):
        self.__rcm_fname = rcm_fname
        self.__clusters_dict = dict()
        self.__read_cluster_map = dict()
        self.__ExtractFromRCM()
        print "File " + self.__rcm_fname + " was successfully parsed"
        print "RCM file consists of " + str(len(self.__clusters_dict)) + " clusters & " + \
              str(len(self.__read_cluster_map)) + " names"

    def __iter__(self):
        for read_name in self.__read_cluster_map:
            yield read_name

    def ContainsRead(self, read_name):
        return read_name in self.__read_cluster_map

    def GetReadCluster(self, read_name):
        return self.__read_cluster_map[read_name]

    def ClustersNumber(self):
        return len(self.__clusters_dict)

    def GetClusterById(self, cluster_id):
        return self.__clusters_dict[cluster_id]

    def GetClusterIds(self):
        return self.__clusters_dict.keys()

def CreateOutputDir(output_dir_name):
    if os.path.exists(output_dir_name):
        shutil.rmtree(output_dir_name)
    os.makedirs(output_dir_name)

def FindReadsExcludedFromNewRepertoire(old_repertoire, new_repertoire, excluded_reads):
    num_excluded = 0
    fhandler = open(excluded_reads, "w")
    for read_name in old_repertoire:
        if not new_repertoire.ContainsRead(read_name):
            num_excluded += 1
            fhandler.write(read_name + "\n")
    fhandler.close()
    print str(num_excluded) + " reads were excluded from the new repertoire"
    print "Excluded reads were written to " + excluded_reads

def FindDisjointPairs(old_repertoire, new_repertoire, disjoint_read_pairs):
    old_cluster_ids = old_repertoire.GetClusterIds()
    fhandler = open(disjoint_read_pairs, "w")
    num_disjoint = 0
    num_all = 0
    for cluster_id in old_cluster_ids:
        cluster = old_repertoire.GetClusterById(cluster_id)
        for j in range(0, len(cluster)):
            for k in range(j + 1, len(cluster)):
                read1 = cluster[j]
                read2 = cluster[k]
                if not new_repertoire.ContainsRead(read1) or not new_repertoire.ContainsRead(read2):
                    continue
                cluster1 = new_repertoire.GetReadCluster(read1)
                cluster2 = new_repertoire.GetReadCluster(read2)
                if cluster1 != cluster2:
                    fhandler.write(read1 + "\t" + read2 + "\t" + str(cluster1) + "\t" + str(cluster2) + "\n")
                    num_disjoint += 1
                num_all += 1
    print str(num_disjoint) + " of " + str(num_all) + \
          " pairs of reads from the same cluster in the old repertoire are " \
          "located in the different clusters of the new repertoire"
    print "Such pairs of reads were written to " + disjoint_read_pairs
    fhandler.close()

def CreateCompressedReadMap(compressed_reads):
    fhandler = open(compressed_reads, "r")
    read_vertex_id_map = dict()
    current_index = 0
    for line in fhandler:
        if line[0] == '>':
            read_name = line.strip()[1:]
            name_end = read_name.find('abundance')
            if name_end == -1:
                print "ERROR: incorrect header " + line + " of file " + compressed_reads
                sys.exit(1)
            read_name = read_name[:(name_end - 1)]
            read_vertex_id_map[read_name] = current_index
            current_index += 1
    fhandler.close()
    return read_vertex_id_map

def FindDisconnectedReadPairs(old_repertoire, read_vertex_id_map, sw_graph, disconnected_pairs):
    fhandler = open(disconnected_pairs, "w")
    old_cluster_ids = old_repertoire.GetClusterIds()
    num_disconnected = 0
    num_all = 0
    for cluster_id in old_cluster_ids:
        cluster = old_repertoire.GetClusterById(cluster_id)
        for j in range(0, len(cluster)):
            for k in range(j + 1, len(cluster)):
                read1 = cluster[j]
                read2 = cluster[k]
                if not read1 in read_vertex_id_map or not read2 in read_vertex_id_map:
                    continue
                index1 = read_vertex_id_map[read1]
                index2 = read_vertex_id_map[read2]
                if not index2 in sw_graph.GetAdjacencyList(index1):
                    fhandler.write(read1 + "\t" + read2 + "\n")
                    num_disconnected += 1
                num_all += 1
    print str(num_disconnected) + " of " + str(num_all) + " pairs of reads from the same cluster in the old repertoire " \
                                                        "were not connected by edge in new SW graph"
    print "Such pairs were written to " + disconnected_pairs
    fhandler.close()

def main():
    if len(sys.argv) != 6:
        print "ERROR: Incorrect input parameters"
        print "python compare_rcm.py old.rcm new.rcm new.graph compressed_reads.fa output_dir"
        sys.exit(1)
    output_dir = sys.argv[5]
    print "Output dir: " + output_dir
    CreateOutputDir(output_dir)
    old_rcm_fname = sys.argv[1]
    new_rcm_fname = sys.argv[2]
    old_repertoire = Repertoire(old_rcm_fname)
    new_repertoire = Repertoire(new_rcm_fname)
    print "===="
    excluded_reads = os.path.join(output_dir, "excluded_reads.txt")
    FindReadsExcludedFromNewRepertoire(old_repertoire, new_repertoire, excluded_reads)
    print "===="
    disjoint_read_pairs = os.path.join(output_dir, "disjoint_read_pairs.txt")
    FindDisjointPairs(old_repertoire, new_repertoire, disjoint_read_pairs)
    print "===="
    graph_fname = sys.argv[3]
    sw_graph = graph_io.MetisGraphReader(graph_fname).ExtractGraph()
    compressed_reads = sys.argv[4]
    read_vertex_id_map = CreateCompressedReadMap(compressed_reads)
    disconnected_pairs = os.path.join(output_dir, "disconnected_pairs.txt")
    FindDisconnectedReadPairs(old_repertoire, read_vertex_id_map, sw_graph, disconnected_pairs)

if __name__ == '__main__':
    main()
