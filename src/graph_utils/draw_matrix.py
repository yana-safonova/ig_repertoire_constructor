import sys
import os
import random
import shutil

import graphics
from graphics import *

import graph_structs
from graph_structs import Graph
from graph_structs import Edge
from graph_structs import Permutation
from graph_structs import Decomposition

import graph_io
from graph_io import MetisGraphReader

class DrawingConfig:
    def __init__(self, graph = Graph(0),
                 permutation = Permutation(0),
                 decomposition = Decomposition(),
                 mode = "",
                 output_img = "",
                 image_size = 600,
                 background_color = "white"):
        self.graph = graph
        self.permutation = permutation
        self.decomposition = decomposition
        self.__mode = mode
        self.output_img = output_img
        self.image_size = image_size
        self.background_color = background_color

    def __CommandLineArgumentsAreCorrect(self, command_arguments):
        if len(command_arguments) == 5 and (command_arguments[2] == 'DEC' or command_arguments[2] == 'PERM'):
            return True
        elif len(command_arguments) == 4 and command_arguments[2] == 'TRIV':
            return True
        return False

    def ModeIsPermutated(self):
        return self.__mode == 'PERM'

    def ModeIsDecomposed(self):
        return self.__mode == 'DEC'

    def ModeIsTrivial(self):
        return self.__mode == 'TRIV'

    def __ModeIsCorrect(self, mode_str):
        return not self.ModeIsPermutated() and \
               not self.ModeIsDecomposed() and \
               not self.ModeIsTrivial()

    def CreateConfigFromCommandLine(self, command_arguments):
        if not self.__CommandLineArgumentsAreCorrect(command_arguments):
            print "Invalid number of input arguments: "
            print "argv[1] - GRAPH filename, " \
                  "argv[2] - mode (PERM / DEC / TRIV), " \
                  "argv[3] - output directory, " \
                  "argv[4] - permutation or decomposition filename for PERM/DEC mode"
            sys.exit(1)
        graph_filename = command_arguments[1]
        self.graph = MetisGraphReader(graph_filename).ExtractGraph()
        mode_string = command_arguments[2]
        if not self.__ModeIsCorrect(mode_string):
            print "Mode string is not correct. Please choose one the following values of mode:"
            print " PERM\t\tdrawing permutated graph,"
            print " DEC\t\tdrawing decomposed graph,"
            print " TRIV\t\tgraph drawing"
            sys.exit(1)
        self.__mode = mode_string
        self.output_img = command_arguments[3]
        if self.ModeIsPermutated():
            self.permutation = Permutation(self.graph.VertexNumber())
            self.permutation.ExtractFromFile(command_arguments[4])
        elif self.ModeIsDecomposed():
            self.decomposition = Decomposition()
            self.decomposition.ExtractFromFile(command_arguments[4])
            self.permutation = Permutation(self.graph.VertexNumber())
            self.permutation.CreatePermutationByDecomposition(self.decomposition)
        elif self.ModeIsTrivial():
            self.permutation = Permutation(self.graph.VertexNumber())
            self.permutation.CreateTrivialPermutation()

    def __str__(self):
        string = "Drawing config:\n"
        string += " - graph with " + str(self.graph.VertexNumber()) + " vertices & " + \
                 str(self.graph.EdgeNumber()) + " edges\n"
        if self.ModeIsPermutated():
            string += " - permutation\n"
        elif self.ModeIsDecomposed():
            string += " - decomposition\n"
        string += " - Output filename: " + self.output_img
        return string

# --------------------- drawing utils ------------------------------

shm_color_dict = {0: color_rgb(41, 209, 37), 1: color_rgb(241, 143, 30), 2: color_rgb(239, 122, 26),
                  3: color_rgb(237, 102, 22), 4: color_rgb(235, 81, 19), 5: color_rgb(234, 61, 15),
                  6: color_rgb(232, 40, 12), 7: color_rgb(230, 20, 8), 8: color_rgb(229, 0, 5)}

weight_color_dict = {1: color_rgb(91, 255, 56), 2: color_rgb(65, 232, 46), 3: color_rgb(41, 209, 37),
                     4: color_rgb(29, 186, 37), 5: color_rgb(22, 163, 40), 6: color_rgb(16, 140, 42),
                     7: color_rgb(11, 117, 41), 8: color_rgb(7, 95, 39)}

dec_color_dict = {0: color_rgb(255, 0, 0), 1: color_rgb(255, 140, 0), 2: color_rgb(255, 255, 0),
                  3: color_rgb(0, 238, 118), 4: color_rgb(0, 191, 255), 5: color_rgb(0, 0, 255),
                  6: color_rgb(238, 0, 238)}

#-----------------------------------------------------------

class GraphDrawer:
    def __init__(self, drawing_config = DrawingConfig()):
        self.__drawing_config = drawing_config
        self.__image = Image(Point(0, 0),
                             self.__drawing_config.image_size,
                             self.__drawing_config.image_size)
        self.__decomposition_color = dict()

    def __PrepareOutputDir(self):
        if os.path.exists(self.__drawing_config.output_dir):
            shutil.rmtree(self.__drawing_config.output_dir)
        os.makedirs(self.__drawing_config.output_dir)

    def __AssignColorForDecomposition(self):
        index = 6
        for i in range(0, self.__drawing_config.decomposition.ClassNumber()):
            self.__decomposition_color[i] = dec_color_dict[i % len(dec_color_dict)]

    def __FillImage(self):
        for i in range(0, self.__drawing_config.image_size):
            for j in range(0, self.__drawing_config.image_size):
                self.__image.setPixel(i, j, self.__drawing_config.background_color)

    def __DrawRectangleOnImage(self, x ,y, x_size, y_size, color):
        for i in range(x, x + x_size):
            for j in range(y, y + y_size):
                self.__image.setPixel(i, j, color)

    def __GetDecompositionColor(self, class_id):
        return self.__decomposition_color[class_id]

    def __GetColorForEdge(self, edge):
        if self.__drawing_config.ModeIsPermutated():
            return weight_color_dict[int(edge.w) % len(weight_color_dict)]
        elif self.__drawing_config.ModeIsDecomposed():
            if self.__drawing_config.decomposition.VerticesFromTheSameClass(edge.v1, edge.v2):
                return self.__GetDecompositionColor(self.__drawing_config.decomposition.GetClassByVertex(edge.v1))
            return "grey"
        return "black"

    def __DrawGraphLessSize(self):
        graph = self.__drawing_config.graph
        rect_size = self.__drawing_config.image_size / graph.VertexNumber() + 1
        for edge in graph:
            if edge.v1 < edge.v2:
                x = int(float(self.__drawing_config.permutation.direct[edge.v1]) /
                        graph.VertexNumber() * self.__drawing_config.image_size)
                y = int(float(self.__drawing_config.permutation.direct[edge.v2]) /
                        graph.VertexNumber() * self.__drawing_config.image_size)
                self.__DrawRectangleOnImage(x, y, rect_size, rect_size, self.__GetColorForEdge(edge))
                self.__DrawRectangleOnImage(y, x, rect_size, rect_size, self.__GetColorForEdge(edge))
        if self.__drawing_config.ModeIsDecomposed():
            for i in range(0, graph.VertexNumber()):
                x = int(float(i) / graph.VertexNumber() * self.__drawing_config.image_size)
                self.__DrawRectangleOnImage(x, x, rect_size, rect_size,
                                            self.__GetColorForEdge(Edge(self.__drawing_config.permutation.reverse[i],
                                                                  self.__drawing_config.permutation.reverse[i],0)))

    def __DrawGraphGreaterSize(self):
        graph = self.__drawing_config.graph
        for edge in graph:
            if not edge.v1 < edge.v2:
                continue
            x = int(float(self.__drawing_config.permutation.direct[edge.v1]) / graph.VertexNumber() *
                    self.__drawing_config.image_size)
            y = int(float(self.__drawing_config.permutation.direct[edge.v2]) / graph.VertexNumber() *
                    self.__drawing_config.image_size)
            self.__image.setPixel(x, y, self.__GetColorForEdge(edge))
            self.__image.setPixel(y, x, self.__GetColorForEdge(edge))

    def __DrawGraph(self):
        if self.__drawing_config.ModeIsDecomposed():
            print "Assignment of colors for decomposition classes..."
            self.__AssignColorForDecomposition()
        print "Filling image background..."
        self.__FillImage()
        print "Drawing graph..."
        if self.__drawing_config.graph.VertexNumber() <= self.__drawing_config.image_size:
            print "Graph size (" + str(self.__drawing_config.graph.VertexNumber()) + ") <= image size (" + \
                  str(self.__drawing_config.image_size) + ")"
            self.__DrawGraphLessSize()
        else:
            print "Graph size (" + str(self.__drawing_config.graph.VertexNumber()) + ") > image size (" + \
                  str(self.__drawing_config.image_size) + ")"
            self.__DrawGraphGreaterSize()
        image_fname = self.__drawing_config.output_img #os.path.join(self.__drawing_config.output_dir, "image.ppm")
        self.__image.save(image_fname)
        print "Matrix image was written to " + image_fname

    def Draw(self):
        #self.__PrepareOutputDir()
        self.__DrawGraph()

#-----------------------------------------------------------

def main():
    drawing_config = DrawingConfig()
    drawing_config.CreateConfigFromCommandLine(sys.argv)
    print drawing_config
    graph_drawer = GraphDrawer(drawing_config)
    graph_drawer.Draw()
    return 0

if __name__ == "__main__":
    main()
