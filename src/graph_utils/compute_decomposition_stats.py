import os
import sys
import draw_matrix

graph_fname = sys.argv[1]
decomposition_fname = sys.argv[2]

graph = draw_matrix.Graph()
graph.ExtractFromGraphFile(graph_fname)

decomposition = draw_matrix.Decomposition()
decomposition.ExtractFromFile(decomposition_fname)

num_classes = decomposition.num_classes
class_fillin = dict()
class_sizes = dict()
for i in range(0, len(decomposition.vertex_color)):
    vertex = i
    color = decomposition.vertex_color[i]
    if color not in class_sizes:
        class_sizes[color] = 0
    class_sizes[color] += 1

for edge in graph.edges:
    if decomposition.vertex_color[edge.v1] != decomposition.vertex_color[edge.v2]:
        continue
    class_id = decomposition.vertex_color[edge.v1]
    if class_id not in class_fillin:
        class_fillin[class_id] = 0
    class_fillin[class_id] += 1

num_singletons = 0
for class_id in class_sizes:
    size = class_sizes[class_id]
    if size != 1:
        fillin = float(class_fillin[class_id]) / size / (size - 1) * 2
        print("Size: " + str(size) + ", edge fill-in: " + str(fillin))
    else:
        num_singletons += 1

print "Num singletons: " + str(num_singletons)



