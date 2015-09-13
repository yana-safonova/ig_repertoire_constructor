__author__ = 'yana'

import os
import sys

# --------------------- Edge ----------------------
class Edge:
    def __init__(self, v1, v2, w = 0):
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

    def __InitializeEdgeList(self):
        for v1 in range(0, len(self.__adjacent_lists)):
            for v2 in self.__adjacent_lists[v1]:
                self.__edges.append(Edge(v1, v2))
        #print str(len(self.__edges)) + " edges were created"
        #sys.exit(1)

    def __init__(self, num_vertices):
        self.__edges = list()
        self.__num_vertices = num_vertices
        self.__num_edges = 0
        self.__InitializeAdjacencyList()

    def AddEdge(self, edge):
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
        if len(self.__edges) == 0:
            self.__InitializeEdgeList()
        for edge in self.__edges:
            yield edge

# --------------------- Decomposition --------------
class Decomposition:
    def __init__(self):
        self.__vertex_color = list()
        self.__num_classes = 0
        self.__decomposition_classes = dict()

    def ExtractFromFile(self, filename):
        fhandler = open(filename, 'r')
        lines = fhandler.readlines()
        for i in range(0, len(lines)):
            set_id = int(lines[i].strip())
            self.__vertex_color.append(set_id)
            self.__num_classes = max(self.__num_classes, set_id)
            if not set_id in self.__decomposition_classes:
                self.__decomposition_classes[set_id] = list()
            self.__decomposition_classes[set_id].append(i)
        self.__num_classes = len(self.__decomposition_classes)
        print "Decomposition into " + str(self.__num_classes) + " classes was extracted from " + filename

    def ClassNumber(self):
        return self.__num_classes

    def __iter__(self):
        for class_id in self.__vertex_color:
            yield class_id

    def VerticesFromTheSameClass(self, v1, v2):
        return self.__vertex_color[v1] == self.__vertex_color[v2]

    def GetClassByVertex(self, vertex):
        return self.__vertex_color[vertex]

    def GetDecompositionClassById(self, class_id):
        return self.__decomposition_classes[class_id]

# --------------------- Permutation ----------------
class Permutation:
    def __init__(self, num_elems):
        self.__num_elems = num_elems
        self.direct = [0] * num_elems # direct[i] shows new index of old vertex i
        self.reverse = [0] * num_elems # reverse[i] shows old index of new vertex i

    def ExtractFromFile(self, filename):
        fhandler = open(filename, "r")
        lines = fhandler.readlines()
        for i in range(0, len(lines)):
            num1 = i
            num2 = int(lines[i].strip())
            self.direct[num1] = num2
            self.reverse[num2] = num1
        fhandler.close()

    def CreateTrivialPermutation(self):
        for i in range(0, self.__num_elems):
            self.direct[i] = i
            self.reverse[i] = i

    def CreatePermutationByDecomposition(self, decompositon = Decomposition()):
        new_index_counter = 0
        for i in range(0, decompositon.ClassNumber()):
            cur_decomposition_class = decompositon.GetDecompositionClassById(i)
            for elem in cur_decomposition_class:
                self.direct[elem] = new_index_counter
                self.reverse[new_index_counter] = elem
                new_index_counter += 1