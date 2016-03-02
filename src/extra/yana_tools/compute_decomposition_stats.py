import sys
import os
import random
import shutil

sys.path.append('../graph_utils_py/')

import graph_structs
from graph_structs import Graph
from graph_structs import Edge
from graph_structs import Permutation
from graph_structs import Decomposition

import graph_io
from graph_io import MetisGraphReader

def main():
    if len(sys.argv) != 3:
        print "ERROR: invalid command line arguments"
        print "python compute_decomposition_stats.py g.graph decomposition.txt"
        sys.exit(1)

    graph_filename = sys.argv[1]
    graph = MetisGraphReader(graph_filename).ExtractGraph()

    dec_filename = sys.argv[2]
    decomposition = Decomposition()
    decomposition.ExtractFromFile(dec_filename)

    fill_dict = dict()
    size_dict = dict()
    for edge in graph:
        if edge.v1 < edge.v2:
            if decomposition.GetClassByVertex(edge.v1) == decomposition.GetClassByVertex(edge.v2):
                class_id = decomposition.GetClassByVertex(edge.v1)
                if not class_id in fill_dict:
                    fill_dict[class_id] = 0
                fill_dict[class_id] += 1

    for i in range(0, graph.VertexNumber()):
        class_id = decomposition.GetClassByVertex(i)
        if not class_id in size_dict:
            size_dict[class_id] = 0
        size_dict[class_id] += 1

    max_class_size = 0
    max_fillin = 0
    avg_fillin = 0
    num_trivial_classes = 0
    num_nt_classes = 0
    for class_id in fill_dict:
        if size_dict[class_id] < 5:
            num_trivial_classes += 1
            continue
        fillin = float(fill_dict[class_id]) / float(size_dict[class_id]) / float(size_dict[class_id] + 1) * 2
        #print str(fillin) + "\t" + str(size_dict[class_id])
        avg_fillin += fillin
        if fillin > max_fillin:
            max_fillin = fillin
        if size_dict[class_id] > max_class_size:
            max_class_size = size_dict[class_id]
        num_nt_classes += 1
    avg_fillin /= float(num_nt_classes)

    print "# decomposition classes:\t\t" + str(len(size_dict))
    print "# trivial decomposition classes:\t" + str(num_trivial_classes)
    print "Size of max class:\t\t\t" + str(max_class_size)
    print "Avg edge fill-in:\t\t\t" + str(avg_fillin)
    print "Max edge fill-in:\t\t\t" + str(max_fillin)

if __name__ == "__main__":
    main()

