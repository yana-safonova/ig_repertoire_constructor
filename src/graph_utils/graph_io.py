__author__ = 'yana'

import os
import sys

import graph_structs

from graph_structs import Graph
from graph_structs import Edge

class MetisGraphReader:
    def __init__(self, filename):
        self.__filename = filename
        if not os.path.exists(self.__filename):
            print "GRAPH file " + self.__filename + " was not found"
            sys.exit(1)
        self.__fhandler = open(self.__filename, 'r')

    def __ParseNumVertices(self, first_line_splits):
        return int(first_line_splits[0])

    def __ParseNumEdges(self, first_line_splits):
        return int(first_line_splits[1])

    def __CreateWeightedGraph(self, filelines):
        for i in range(1, len(filelines)):
            splits = filelines[i].strip().split()
            if len(splits) % 2 != 0:
                print "Incorrect line of GRAPH file " + self.__filename + " contains odd number of elements:"
                print filelines[i]
                sys.exit(1)
            for j in range(0, len(splits) / 2):
                neigh = int(splits[j * 2]) - 1
                weight = float(splits[j * 2 + 1])
                self.__graph.AddEdge(Edge(i - 1, neigh, weight))

    def __CreateUnweightedGraph(self, filelines):
        for i in range(1, len(filelines)):
            splits = filelines[i].strip().split()
            for j in range(0, len(splits)):
                neigh = int(splits[j]) - 1
                weight = 1.0
                self.__graph.AddEdge(Edge(i - 1, neigh, weight))

    def __GraphIsUnweighted(self, first_line_splits):
        return len(first_line_splits) == 2

    def __GraphIsWeighted(self, first_line_splits):
        if len(first_line_splits) < 3:
            return False
        return first_line_splits[2] == '001'

    def ExtractGraph(self):
        lines = self.__fhandler.readlines()
        first_line_splits = lines[0].strip().split()
        self.__graph = Graph(self.__ParseNumVertices(first_line_splits))
        if self.__GraphIsUnweighted(first_line_splits):
            #print "Unweighted graph reader was chosen"
            self.__CreateUnweightedGraph(lines)
        elif self.__GraphIsWeighted(first_line_splits):
            #print "Weighted graph reader was chosen"
            self.__CreateWeightedGraph(lines)
        print "Graph with " + str(self.__graph.VertexNumber()) + " & " + str(self.__graph.EdgeNumber()) + \
              " edges was extracted from " + self.__filename
        self.__fhandler.close()
        return self.__graph

###########################

class MetisGraphWriter:
    def __init__(self, graph):
        self.__graph = graph

    def WriteInFile(self, filename):
        fhandler = open(filename, "w")
        fhandler.write(str(self.__graph.VertexNumber()) + "\t" + str(self.__graph.EdgeNumber()) + "\n")
        for v1 in range(0, self.__graph.VertexNumber()):
            if not self.__graph.VertexIsIsolated(v1):
                adj_list = self.__graph.GetAdjacencyList(v1)
                for v2 in adj_list:
                    fhandler.write(str(v2 + 1) + " ")
            fhandler.write("\n")
        print "Graph with " + str(self.__graph.VertexNumber()) + " vertices & " + str(self.__graph.EdgeNumber()) + \
              " edges was written to " + filename
        fhandler.close()

###################################################

class MtxGraphReader:
    def __init__(self, filename):
        self.__filename = filename
        if not os.path.exists(self.__filename):
            print "MTX file " + self.__filename + " was not found"
            sys.exit(1)
        self.__fhandler = open(self.__filename, 'r')

    def __GetIndexOfHeaderLine(self, filelines):
        for i in range(0, len(filelines)):
            if filelines[i][0] != '%':
                return i
        return -1

    def __MatrixIsSymmetric(self, header_line_splits):
        return header_line_splits[0] == header_line_splits[1]

    def __ParseNumVertices(self, header_line_splits):
        return int(header_line_splits[1])

    def __ParseNumEdges(self, header_line_splits):
        return int(header_line_splits[2])

    def ExtractGraph(self):
        lines = self.__fhandler.readlines()
        header_index = self.__GetIndexOfHeaderLine(lines)
        header_line_splits = lines[header_index].strip().split()
        if not self.__MatrixIsSymmetric(header_line_splits):
            print "WARN: graph " + self.__filename + " is not symmetric. Undirected graph can not be constructed"
            sys.exit(1)
        graph = Graph(self.__ParseNumVertices(header_line_splits))
        for i in range(header_index + 1, len(lines)):
            splits = lines[i].strip().split()
            if splits[0] == splits[1]:
                continue
            weight = 1
            if len(splits) > 2:
                weight = float(splits[2])
            v1 = int(splits[0]) - 1
            v2 = int(splits[1]) - 1
            graph.AddEdge(Edge(v1, v2, weight))
        print "Graph with " + str(graph.VertexNumber()) + " vertices & " + str(graph.EdgeNumber()) + \
              " edges was extracted from " + self.__filename
        self.__fhandler.close()
        return graph

###########################

class MtxGraphWriter:
    def __init__(self, graph):
        self.__graph = graph

    def WriteInFile(self, filename):
        fhandler = open(filename, 'w')
        fhandler.write("%%MatrixMarket matrix coordinate real general\n")
        fhandler.write("% Created in laboratory \"Center for Algorithmic Biotechnology\", "
                       "Saint-Petersburg State Unversity, Russia\n")
        fhandler.write("% author: Yana Safonova\n")
        fhandler.write(str(self.__graph.VertexNumber()) + "\t" + str(self.__graph.EdgeNumber()) + "\n")
        for edge in self.__graph:
            fhandler.write(str(edge.v1) + "\t" + str(edge.v2) + '\t' + str(edge.w) + '\n')
        fhandler.close()