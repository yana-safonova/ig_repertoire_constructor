import os
import sys

input_file = sys.argv[1]
output_file = sys.argv[1] + ".graph"

class Edge:
    vertex = 0
    weight = 0

    def __init__(self, vertex, weight):
        self.vertex = vertex
        self.weight = weight

class InputGraph:
    adj_list = dict()
    num_vertices = 0
    num_edges = 0
    w_dict = {0:4, 1:3, 2:2, 3:1}

    def __init__(self):
        self.adj_list = dict()
        self.num_vertices = 0
        self.num_edges = 0

    def AddEdge(self, v1, v2, w):
        edge1 = Edge(v1, w)
        edge2 = Edge(v2, w)
        if v2 not in self.adj_list:
            self.adj_list[v2] = list()
        self.adj_list[v2].append(edge1)
        if v1 not in self.adj_list:
            self.adj_list[v1] = list()
        self.adj_list[v1].append(edge2)

    def ExtractFromFile(self, filename):
        fhandler = open(filename, "r")
        lines = fhandler.readlines()
        self.num_vertices = int(lines[0].strip()[2:])
        self.num_edges = len(lines) - 1
        for i in range(1, len(lines)):
            splits = lines[i].strip().split()
            v1 = int(splits[0]) + 1
            v2 = int(splits[1]) + 1
            w = int(splits[2]) #self.w_dict[int(splits[2])]
            self.AddEdge(v1, v2, w)
        fhandler.close()

def WriteOutputWeightedGraph(input_graph, output_file):
    fhandler = open(output_file, "w")
    fhandler.write(str(input_graph.num_vertices) + " " + str(input_graph.num_edges) + " 001\n")
    for i in input_graph.adj_list:
        adj_list = input_graph.adj_list[i]
        for p in adj_list:
            fhandler.write(str(p.vertex) + " " + str(p.weight) + " ")
        fhandler.write("\n")

def WriteOutputSimpleGraph(input_graph, output_file):
    fhandler = open(output_file, "w")
    fhandler.write(str(input_graph.num_vertices) + " " + str(input_graph.num_edges) + "\n")
    for i in input_graph.adj_list:
        adj_list = input_graph.adj_list[i]
        for p in adj_list:
            fhandler.write(str(p.vertex) + " ")
        fhandler.write("\n")

input_graph = InputGraph()
input_graph.ExtractFromFile(input_file)

if input_graph.num_vertices < 100:
    print "Error: Hamming graph is too small"
    sys.exit(1)
else:
    WriteOutputWeightedGraph(input_graph, output_file)
    #WriteOutputSimpleGraph(input_graph, output_file)
    print "Hamming graph was written to " + output_file
