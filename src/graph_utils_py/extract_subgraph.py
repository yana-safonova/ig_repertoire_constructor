import sys
import os
import draw_matrix

input_graph = sys.argv[1]
tau_threshold = int(sys.argv[2])

graph = draw_matrix.Graph()
graph.ExtractFromGraphFile(input_graph)

output_graph = input_graph + "_tau_" + str(tau_threshold)
output_fhandler = open(output_graph, "w")
output_fhandler.write("n=" + str(graph.num_vertices)  + "\n")
for edge in graph.edges:
    if edge.w <= tau_threshold:
        output_fhandler.write(str(edge.v1) + " " + str(edge.v2) + " " + str(edge.w) + " " + str(edge.shm) + "\n")
output_fhandler.close()
