import os
import sys
import shutil

graph_utils = "src/graph_utils"
sys.path.append(graph_utils)

import graph_structs
from graph_structs import Graph
import graph_io

from Bio import SeqIO

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

################################################################

class NewRepertoire:
    def __CreateCompressedMap(self, compressed_fa_fname):
        fhandler = open(compressed_fa_fname, "r")
        current_index = 0
        for line in fhandler:
            if line[0] == '>':
                read_name = line.strip()[1:]
                name_end = read_name.find('abundance')
                if name_end == -1:
                    print "ERROR: incorrect header " + line + " of file " + compressed_fa_fname
                    sys.exit(1)
                read_name = read_name[:(name_end - 1)]
                self.readname_vertexid[read_name] = current_index
                self.vertexid_readname.append(read_name)
                current_index += 1
        fhandler.close()

    def __ReadCropppedIds(self, cropped_fa):
        for record in SeqIO.parse(open(cropped_fa), 'fasta'):
            self.cropped_read_names.append(record.id)
        print str(len(self.cropped_read_names)) + "  names of cropped reads were extracted from " + cropped_fa

    def __ReadCroppedReadMap(self, read_map):
        fhandler = open(read_map, "r")
        compressed_indices = list()
        for line in fhandler:
            compressed_indices.append(int(line.strip()))
        fhandler.close()
        if len(self.cropped_read_names) != len(compressed_indices):
            print "ERROR: Cropped read map & cropped FA are not consistent"
            sys.exit(1)
        for i in range(0, len(self.cropped_read_names)):
            read_cropped = self.cropped_read_names[i]
            new_index = compressed_indices[i]
            read_compressed = self.vertexid_readname[new_index]
            self.cropped_read_map[read_cropped] = read_compressed

    def __init__(self, compressed_fa_fname, graph_fname, read_map, croppped_fa):
        self.sw_graph = graph_io.MetisGraphReader(graph_fname).ExtractGraph()
        self.readname_vertexid = dict()
        self.vertexid_readname = list()
        self.__CreateCompressedMap(compressed_fa_fname)
        self.cropped_read_names = list()
        self.__ReadCropppedIds(croppped_fa)
        self.cropped_read_map = dict()
        self.__ReadCroppedReadMap(read_map)

################################################################

class OldRepertoire:
    def __CreateReadVertexMaps(self, old_decomposition):
        fhandler = open(old_decomposition, "r")
        current_index = 0
        for line in fhandler:
            splits = line.strip().split()
            read_name = splits[1]
            self.vertexid_readname.append(read_name)
            self.readname_vertexid[read_name] = current_index
            current_index += 1
        fhandler.close()

    def __CreateReadMap(self, old_read_map):
        fhandler = open(old_read_map, "r")
        for line in fhandler:
            splits = line.strip().split()
            self.read_map[splits[0]] = splits[1]
        fhandler.close()

    def __init__(self, old_graph, old_decomposition):
        self.sw_graph = graph_io.MetisGraphReader(old_graph).ExtractGraph()
        self.readname_vertexid = dict()
        self.vertexid_readname = list()
        self.__CreateReadVertexMaps(old_decomposition)
        #self.read_map = dict()
        #self.__CreateReadMap(old_read_map)

################################################################

def CreateStructFromInputDir(new_igrepcon_dir):
    if not os.path.exists(new_igrepcon_dir):
        print "ERROR: Directory " + new_igrepcon_dir + " was not found"
        sys.exit(1)
    compressed_fa = os.path.join(new_igrepcon_dir, "compressed.fa")
    if not os.path.exists(compressed_fa):
        print "ERROR: Compressed reads " + compressed_fa + " was not found"
        sys.exit(1)
    sw_graph = os.path.join(new_igrepcon_dir, "sw.graph")
    if not os.path.exists(sw_graph):
        print "ERROR: SW graph " + sw_graph + " was not found"
        sys.exit(1)
    read_map = os.path.join(new_igrepcon_dir, "map.txt")
    if not os.path.exists(read_map):
        print "ERROR: Read map " + read_map + " was not found"
        sys.exit(1)
    cropped_fa = os.path.join(new_igrepcon_dir, "vj_finder/cropped.fa")
    if not os.path.exists(cropped_fa):
        print "ERROR: Cropped reads " + cropped_fa + " was not found"
        sys.exit(1)
    return NewRepertoire(compressed_fa, sw_graph, read_map, cropped_fa)

################################################################

def CompareGraph(old_repertoire, new_repertoire):
    old_graph = old_repertoire.sw_graph
    num_disconnected = 0
    num_all = 0
    for i in range(0, old_graph.VertexNumber()):
        adj_list = old_graph.GetAdjacencyList(i)
        for j in adj_list:
            if j > i:
                print "=="
                old_read1 = old_repertoire.vertexid_readname[i]
                old_read2 = old_repertoire.vertexid_readname[j]
                print "Old reads: " + old_read1 + " " + old_read2
                if not old_read1 in new_repertoire.cropped_read_map or not old_read2 in new_repertoire.cropped_read_map:
                    print "Read " + old_read1 + " or " + old_read2 + " was discarded"
                    continue
                read1 = new_repertoire.cropped_read_map[old_read1]
                read2 = new_repertoire.cropped_read_map[old_read2]
                print "New reads: " + read1 + " " + read2
                if not read1 in new_repertoire.readname_vertexid or not read2 in new_repertoire.readname_vertexid:
                    print "ERROR: reads " + read1 + " " + read2 + " were not found"
                    sys.exit(1)
                new_index1 = new_repertoire.readname_vertexid[read1]
                new_index2 = new_repertoire.readname_vertexid[read2]
                print "New indices: " + str(new_index1) + " " + str(new_index2)
                if read1 == read2:
                    print "Reads were compressed together"
                    continue
                if not new_index2 in new_repertoire.sw_graph.GetAdjacencyList(new_index1):
                    print read1 + "\t" + read2 + " were not connected by edge in the new graph"
                    num_disconnected += 1
                num_all += 1
    print str(num_disconnected) + " of " + str(num_all) + " edges were missed in the new graph with respect to the old graph"

def main():
    if len(sys.argv) != 4:
        print "ERROR: Incorrect input parameters"
        print "python compare_rcm.py new_igrepcon_dir old.graph old_decomposition.txt"
        sys.exit(1)
    new_igrepcon_dir = sys.argv[1]
    new_repertoire = CreateStructFromInputDir(new_igrepcon_dir)
    print "===="
    old_graph = sys.argv[2]
    old_decomposition = sys.argv[3]
    old_repertoire = OldRepertoire(old_graph, old_decomposition)
    print "===="
    CompareGraph(old_repertoire, new_repertoire)

if __name__ == '__main__':
    main()
