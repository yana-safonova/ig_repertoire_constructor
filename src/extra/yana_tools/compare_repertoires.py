import os
import sys
import shutil

graph_utils = "src/graph_utils_py"
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

    def __ReadDecomposition(self, decomposition):
        fhandler = open(decomposition, "r")
        for line in fhandler:
            self.decomposition.append(int(line.strip()))
        fhandler.close()
        print "Decomposition of length " + str(len(self.decomposition)) + " was extracted from " + decomposition

    def __CreateReadDecomposition(self):
        for i in range(0, len(self.decomposition)):
            read_name = self.vertexid_readname[i]
            cluster_id = self.decomposition[i]
            self.read_decomposition[read_name] = cluster_id

    def __init__(self, compressed_fa_fname, graph_fname, read_map, croppped_fa, decomposition):
        self.sw_graph = graph_io.MetisGraphReader(graph_fname).ExtractGraph()
        self.readname_vertexid = dict()
        self.vertexid_readname = list()
        self.__CreateCompressedMap(compressed_fa_fname)
        self.cropped_read_names = list()
        self.__ReadCropppedIds(croppped_fa)
        self.cropped_read_map = dict()
        self.__ReadCroppedReadMap(read_map)
        self.decomposition = list()
        self.__ReadDecomposition(decomposition)
        self.read_decomposition = dict()
        self.__CreateReadDecomposition()

################################################################
class OldRepertoireFiles:
    def __init__(self, output_dir):
        self.rcm = os.path.join(output_dir, "constructed_repertoire.rcm")
        if not os.path.exists(self.rcm):
            print "ERROR: RCM file " + self.rcm + " was not found"
            sys.exit(1)
        self.hamming_graph_dir = os.path.join(output_dir, "hamming_graphs_tau_3")
        if not os.path.exists(self.hamming_graph_dir) or not os.path.isdir(self.hamming_graph_dir):
            print "ERROR: Hamming graph directory " + self.hamming_graph_dir + " was not found"
            sys.exit(1)
        print "Hamming graphs directory: " + self.hamming_graph_dir
        self.dense_sgraphs_dir = os.path.join(output_dir, "dense_subgraphs")
        if not os.path.exists(self.dense_sgraphs_dir) or not os.path.isdir(self.dense_sgraphs_dir):
            print "ERROR: Dense sgraph decomposition directory " + self.dense_sgraphs_dir + " was not found"
        print "Dense subgraphs and read maps directory: " + self.dense_sgraphs_dir

class OldGraphDecomposition:
    def __CreateReadVertexMaps(self, decomposition_fname):
        fhandler = open(decomposition_fname, "r")
        current_index = 0
        for line in fhandler:
            splits = line.strip().split()
            if len(splits) != 3:
                print "Line " + line.strip() + " of file " + decomposition_fname + " is incorrect"
                sys.exit(1)
            read_name = splits[1]
            self.vertexid_readname.append(read_name)
            self.readname_vertexid[read_name] = current_index
            current_index += 1
            self.vertex_decomposition[read_name] = int(splits[0])
        print "Dense subgraph decomposition for " + str(len(self.readname_vertexid)) + \
              " vertices was extracted from " + decomposition_fname
        fhandler.close()

    def __CreateReadNameMap(self, read_map_fname):
        fhandler = open(read_map_fname, "r")
        for line in fhandler:
            splits = line.strip().split()
            self.read_name_map[splits[0]] = splits[1]
        print str(len(self.read_name_map)) + " read pairs were extracted from " + read_map_fname
        fhandler.close()

    def __init__(self, graph_fname, decomposition_fname, read_map_fname):
        self.sw_graph = graph_io.MetisGraphReader(graph_fname).ExtractGraph()
        self.readname_vertexid = dict()
        self.vertexid_readname = list()
        self.vertex_decomposition = dict()
        self.__CreateReadVertexMaps(decomposition_fname)
        self.read_name_map = dict()
        self.__CreateReadNameMap(read_map_fname)

class OldRepertoire:
    def __init__(self, output_dir):
        self.repertoire_files = OldRepertoireFiles(output_dir)
        self.repertoire = Repertoire(self.repertoire_files.rcm)
        self.__hamming_graph_fnames = [os.path.join(self.repertoire_files.hamming_graph_dir, f) for f
                                       in os.listdir(self.repertoire_files.hamming_graph_dir)
                                       if os.path.isfile(os.path.join(self.repertoire_files.hamming_graph_dir, f))]
        print str(len(self.__hamming_graph_fnames)) + " Hamming graphs were found in " + \
              self.repertoire_files.hamming_graph_dir

    def __iter__(self):
        suffix_iperm = ".iperm"
        for hg_fname in self.__hamming_graph_fnames:
            if hg_fname[len(hg_fname) - len(suffix_iperm):] == suffix_iperm:
                continue
            yield hg_fname

    def __GetHammingGraphId(self, hg_fname):
        hg_fname = os.path.basename(hg_fname)
        splits = hg_fname.split('.')
        if len(splits) != 2:
            print "ERROR: Hamming graph fname " + hg_fname + " is not correct"
            sys.exit(1)
        splits = splits[0].split("_")
        if len(splits) != 3:
            print "ERROR: Hamming graph fname " + hg_fname + " is not correct"
            sys.exit(1)
        return splits[2]

    def GetDecompositionFname(self, hg_fname):
        graph_id = self.__GetHammingGraphId(hg_fname)
        return os.path.join(self.repertoire_files.dense_sgraphs_dir,
                            "dense_sgraph_decomposition_" + str(graph_id) + ".txt")

    def DecompositionForGraphExists(self, hg_fname):
        return os.path.exists(self.GetDecompositionFname(hg_fname))

    def GetReadMapFname(self, hg_fname):
        graph_id = self.__GetHammingGraphId(hg_fname)
        return os.path.join(self.repertoire_files.dense_sgraphs_dir,
                            "read_map_" + str(graph_id) + ".txt")

    def ReadMapForGraphExists(self, hg_fname):
        return os.path.join(self.GetReadMapFname(hg_fname))

    def CreateSubgraphDecomposition(self, hg_fname, dec_fname, read_map_fname):
        return OldGraphDecomposition(hg_fname, dec_fname, read_map_fname)

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
    decomposition = os.path.join(new_igrepcon_dir, "dense_sgraph_finder/dense_subgraphs.txt")
    if not os.path.exists(decomposition):
        print "ERROR: Dense subgraph decomposition " + decomposition + " was not found"
        sys.exit(1)
    return NewRepertoire(compressed_fa, sw_graph, read_map, cropped_fa, decomposition)

################################################################
class ReadPair:
    def __init__(self, read1, read2):
        self.read1 = read1
        self.read2 = read2

def CompareGraph(old_graph_dec, new_repertoire):
    old_graph = old_graph_dec.sw_graph
    num_all = 0
    disconnected_pairs = list()
    for i in range(0, old_graph.VertexNumber()):
        adj_list = old_graph.GetAdjacencyList(i)
        for j in adj_list:
            if j > i:
                old_read1 = old_graph_dec.vertexid_readname[i]
                old_read2 = old_graph_dec.vertexid_readname[j]
                if not old_read1 in new_repertoire.cropped_read_map or not old_read2 in new_repertoire.cropped_read_map:
                    #print "Read " + old_read1 + " or " + old_read2 + " was discarded"
                    continue
                read1 = new_repertoire.cropped_read_map[old_read1]
                read2 = new_repertoire.cropped_read_map[old_read2]
                #print "New reads: " + read1 + " " + read2
                if not read1 in new_repertoire.readname_vertexid or not read2 in new_repertoire.readname_vertexid:
                    #print "ERROR: reads " + read1 + " " + read2 + " were not found"
                    sys.exit(1)
                new_index1 = new_repertoire.readname_vertexid[read1]
                new_index2 = new_repertoire.readname_vertexid[read2]
                #print "New indices: " + str(new_index1) + " " + str(new_index2)
                if read1 == read2:
                    #print "Reads were compressed together"
                    continue
                if not new_index2 in new_repertoire.sw_graph.GetAdjacencyList(new_index1):
                    print "WARN: " + read1 + "\t" + read2 + " were not connected by edge in the new graph"
                    disconnected_pairs.append(ReadPair(read1, read2))
                num_all += 1
    if len(disconnected_pairs) != 0:
        print "WARN: " + str(len(disconnected_pairs)) + " of " + str(num_all) + \
              " edges were missed in the new graph with respect to the old graph"
    return disconnected_pairs

def CompareGraphs(old_repertoire, new_repertoire, output_fname):
    num_disconnected = 0
    output_fhandler = open(output_fname, "w")
    for hgraph in old_repertoire:
        print "Analysis of Hamming graph " + hgraph + " starts"
        dec = old_repertoire.GetDecompositionFname(hgraph)
        if not old_repertoire.DecompositionForGraphExists(hgraph):
            print "WARN: Dense subgraph decomposition for Hamming graph " + hgraph + " was not found"
            continue
        read_map = old_repertoire.GetReadMapFname(hgraph)
        if not old_repertoire.ReadMapForGraphExists(hgraph):
            print "WARN: Read map for Hamming graph " + hgraph + " was not found"
            continue
        graph_decomposition = old_repertoire.CreateSubgraphDecomposition(hgraph, dec, read_map)
        disconnected_pairs = CompareGraph(graph_decomposition, new_repertoire)
        # output disconnected reads
        for read_pair in disconnected_pairs:
            output_fhandler.write(read_pair.read1 + "\t" + read_pair.read2 + "\n")
        num_disconnected += len(disconnected_pairs)
        print "Analysis of Hamming graph " + hgraph + " ends"
        print "===="
    output_fhandler.close()
    print str(num_disconnected) + " edges missed in new graph compared to new graph were written to " + output_fname

#######################################
def CreateClusterDict(decomposition):
    cluster_dict = dict()
    for read in decomposition:
        cluster = decomposition[read]
        if not cluster in cluster_dict:
            cluster_dict[cluster] = list()
        cluster_dict[cluster].append(read)
    return cluster_dict

def CompareDecomposition(cluster_dict, new_repertoire):
    num_disconnected = 0
    for cluster in cluster_dict:
        read_list = cluster_dict[cluster]
        for i in range(0, len(read_list)):
            for j in range(i + 1, len(read_list)):
                old_read1 = read_list[i]
                old_read2 = read_list[j]
                if not old_read1 in new_repertoire.cropped_read_map or not old_read2 in new_repertoire.cropped_read_map:
                    continue
                new_read1 = new_repertoire.cropped_read_map[old_read1]
                new_read2 = new_repertoire.cropped_read_map[old_read2]
                if new_read1 == new_read2:
                    continue
                new_index1 = new_repertoire.readname_vertexid[new_read1]
                new_index2 = new_repertoire.readname_vertexid[new_read2]
                if not new_index2 in new_repertoire.sw_graph.GetAdjacencyList(new_index1):
                    continue
                if new_repertoire.read_decomposition[new_read1] != new_repertoire.read_decomposition[new_read2]:
                    print "WARN: reads " + new_read1 + " " + new_read2 + " joined in the same cluster of the old graph were disconnected"
                    num_disconnected += 1
    if num_disconnected != 0:
        print str(num_disconnected) + " pairs of reads were disconnected in different dense subgraphs"
    return num_disconnected

def CompareDecompositions(old_repertoire, new_repertoire, output_fname):
    num_disconnected = 0
    for hgraph in old_repertoire:
        print "Analysis of Hamming graph " + hgraph + " starts"
        dec = old_repertoire.GetDecompositionFname(hgraph)
        if not old_repertoire.DecompositionForGraphExists(hgraph):
            print "WARN: Dense subgraph decomposition for Hamming graph " + hgraph + " was not found"
            continue
        read_map = old_repertoire.GetReadMapFname(hgraph)
        if not old_repertoire.ReadMapForGraphExists(hgraph):
            print "WARN: Read map for Hamming graph " + hgraph + " was not found"
            continue
        graph_decomposition = old_repertoire.CreateSubgraphDecomposition(hgraph, dec, read_map)
        cluster_dict = CreateClusterDict(graph_decomposition.vertex_decomposition)
        num_disconnected += CompareDecomposition(cluster_dict, new_repertoire)
        print "Analysis of Hamming graph " + hgraph + " ends"
        print "===="
    print str(num_disconnected) + " reads were disconnected in different dense subgraphs"

#######################################
def PrepareOutputDir(output_dir):
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

#######################################
def main():
    if len(sys.argv) != 4:
        print "ERROR: Incorrect input parameters"
        print "python compare_rcm.py new_igrepcon_dir old_igrepcon_dir output_dir"
        sys.exit(1)
    new_igrepcon_dir = sys.argv[1]
    print "New repertoire constructor directory: " + new_igrepcon_dir
    new_repertoire = CreateStructFromInputDir(new_igrepcon_dir)
    print "===="
    old_igrepcon_dir = sys.argv[2]
    print "Old repertoire constructor directory: " + old_igrepcon_dir
    old_repertoire = OldRepertoire(old_igrepcon_dir)
    print "===="
    output_dir = sys.argv[3]
    print "Output directory: " + output_dir
    PrepareOutputDir(output_dir)
    print "===="
    print "Graph comparison starts"
    output_fname = os.path.join(output_dir, "disconnected_reads.txt")
    CompareGraphs(old_repertoire, new_repertoire, output_fname)
    print "===="
    print "Dense subgraph decomposition starts"
    output_fname2 = os.path.join(output_dir, "disjoint_pairs.txt")
    CompareDecompositions(old_repertoire, new_repertoire, output_fname2)

if __name__ == '__main__':
    main()
