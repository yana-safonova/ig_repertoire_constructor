import os
import sys

# --------------------- Edge ----------------------
class Edge:
    v1 = 0
    v2 = 0
    w = 0
    shm = 0
    
    def __init__(self, v1, v2, w, shm):
        self.v1 = v1
        self.v2 = v2
        self.w = w
        self.shm = shm

# --------------------- Graph ---------------------
class Graph:
    edges = list()
    num_vertices = 0
    name = ""

    def __init(self):
        self.edges = list()

    def ParseNumVertices(self, first_line):
        self.num_vertices = int(first_line.strip()[2:])

    def ExtractFromFile(self, filename):
        self.name = filename
        fhandler = open(filename, "r")
        lines = fhandler.readlines()
        self.ParseNumVertices(lines[0])
        for i in range(1, len(lines)):
            splits = lines[i].strip().split()
            edge = Edge(int(splits[0]), int(splits[1]), int(splits[2]), int(splits[3]))
            self.edges.append(edge)
        fhandler.close()
# ------------------------------------------------------

input_graph = sys.argv[1]
graph = Graph()
graph.ExtractFromFile(input_graph)

disjoint_set = list()
for i in range(0, graph.num_vertices):
    disjoint_set.append(i)

for edge in graph.edges:
    if edge.w == 0:
        collapsed_vertex = max(edge.v1, edge.v2)
        main_vertex = min(edge.v1, edge.v2)
        disjoint_set[collapsed_vertex] = main_vertex

num_main_vertices = 0
new_ids = dict()
for i in range(0, len(disjoint_set)):
    if disjoint_set[i] == i:
        new_ids[i] = num_main_vertices
        num_main_vertices += 1

if num_main_vertices < 100:
    print "Collapsed graph " + input_graph + " is too small"
    sys.exit(0)

output_graph = input_graph + "_nsize_" + str(num_main_vertices)
output_fhandler = open(output_graph, 'w')
output_fhandler.write("n=" + str(num_main_vertices) + "\n")
for edge in graph.edges:
    if disjoint_set[edge.v1] == edge.v1 and disjoint_set[edge.v2] == edge.v2:
        output_fhandler.write(str(new_ids[edge.v1]) + " " + str(new_ids[edge.v2]) + " " + str(edge.w)  + " " + str(edge.shm) + "\n")
output_fhandler.close()

print "Collapsed graph was written from " + input_graph + " to " + output_graph
