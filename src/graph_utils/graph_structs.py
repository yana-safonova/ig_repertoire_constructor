__author__ = 'yana'

import os
import sys

# --------------------- Edge ----------------------
class Edge:
    def __init__(self, v1, v2, w):
        self.v1 = v1
        self.v2 = v2
        self.w = w

# --------------------- Graph ---------------------
class Graph:
    def __InitializeAdjacencyList(self):
        self.__adjacent_lists = list()
        for i in range(0, self.__num_vertices):
            self.__adjacent_lists.append(set())

    def __ComputeEdgeNumber(self):
        for adj_list in self.__adjacent_lists:
            self.__num_edges += len(adj_list)
        self.__num_edges /= 2

    def __init__(self, num_vertices):
        self.__edges = list()
        self.__num_vertices = num_vertices
        self.__num_edges = 0
        self.__InitializeAdjacencyList()

    def AddEdge(self, edge):
        self.__edges.append(edge)
        self.__adjacent_lists[edge.v1].add(edge.v2)
        self.__adjacent_lists[edge.v2].add(edge.v1)

    def GetAdjacencyList(self, vertex):
        return self.__adjacent_lists[vertex]

    def VertexIsIsolated(self, vertex):
        if vertex < self.__num_vertices and len(self.__adjacent_lists[vertex]) == 0:
            return True
        return False

    def VertexNumber(self):
        return self.__num_vertices

    def EdgeNumber(self):
        if self.__num_edges == 0:
            self.__ComputeEdgeNumber()
        return self.__num_edges

    def __iter__(self):
        for edge in self.__edges:
            yield edge

# --------------------- Permutation ----------------
class Permutation:
    def __init__(self, num_elems):
        self.direct = [0] * num_elems
        self.reverse = [0] * num_elems

    def ExtractFromFile(self, filename):
        fhandler = open(filename, "r")
        lines = fhandler.readlines()
        for i in range(0, len(lines)):
            num1 = i
            num2 = int(lines[i].strip())
            self.direct[num1] = num2
            self.reverse[num2] = num1
        fhandler.close()

# --------------------- Decomposition --------------
class Decomposition:
    vertex_color = list()
    class_color = dict()
    num_classes = 0

    def __init__(self):
        self.vertex_color = list()
        self.class_color = dict()
        self.num_classes = 0

    def ExtractFromFile(self, filename):
        fhandler = open(filename, 'r')
        lines = fhandler.readlines()
        for l in lines:
            set_id = int(l.strip())
            self.vertex_color.append(set_id)
            self.class_color[set_id] = 'red'
        self.num_classes = len(self.class_color)
        print "Decomposition into " + str(self.num_classes) + " classes was extracted from " + filename
