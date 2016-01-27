import os
import sys

import graph_structs
from graph_structs import Graph
from graph_structs import Edge
import graph_io

def main():
    if len(sys.argv) < 3:
        print "python mtx_to_graph.py input_graph.mtx output_graph.graph"
        sys.exit(1)
    mtx_fname = sys.argv[1]
    graph_fname = sys.argv[2]
    mtx_graph_reader = graph_io.MtxGraphReader(mtx_fname)
    graph = mtx_graph_reader.ExtractGraph()
    graph_writer = graph_io.MetisGraphWriter(graph)
    graph_writer.WriteInFile(graph_fname)

if __name__ == "__main__":
    main()

