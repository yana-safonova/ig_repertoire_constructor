import os
import sys
import draw_matrix

def CreateDotFile(graph, decomposition, output_graph_fname):
    fhandler = open(output_graph_fname, "w")
    fhandler.write("graph g {" + "\n")
    # writing vertices
    for i in range(0, graph.num_vertices):
        color = decomposition.vertex_color[i]
        fhandler.write(str(i) + ' [fillcolor = ' + str(color) + "];\n")
    # writing edges
    for e in graph.edges:
        fhandler.write(str(e.v1) + ' -- ' + str(e.v2) + ";\n")
    fhandler.write('\n')

input_graph_fname = sys.argv[1]
decomposition_fname = sys.argv[2]
#permutation_fname = sys.argv[3]
output_fname = sys.argv[3]

graph = draw_matrix.Graph()
graph.ExtractFromGraphFile(input_graph_fname)

decomposition = draw_matrix.Decomposition()
decomposition.ExtractFromFile(decomposition_fname)
    
#permutation = draw_matrix.Permutation(graph.num_vertices)
#permutation.ExtractFromFile(permutation_fname)
#
#for i in range(0, graph.num_vertices):
#    old_index = permutation.reverse[i]
#    print "Vertex: " + str(old_index) + ", color: " + str(decomposition.vertex_color[old_index]) 

CreateDotFile(graph, decomposition, output_fname)
