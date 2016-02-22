import os
import sys
import igraph
import ntpath
#os.system("export LC_ALL=en_US.UTF-8")
#os.system("export LANG=en_US.UTF-8")
#from graph_tool.all import *

def get_db_number(vertices_labels):
    for i in range(0, len(vertices_labels)):
        if len(vertices_labels[i]) == 2:
            return i
    return -1

def get_color_by_label(label):
    if label[0] == "H":
            return "#4169E1" # "blue"
    elif label[0] == "K":
        return "#ff3333" #"red"
    elif label[0] == "L":
        return "#228B22" # "green"
    return "white" #"black"

def get_shape_by_label(label):
    if label[0] == "H" or label[0] == "K" or label[0] == "L":
        return "circle"
    return "rectangle"

def get_text_position(label):
    if label[0] == "H" or label[0] == "K" or label[0] == "L":
        return -1
    return -5

def draw_with_igraph(vertices_labels, edges, output_pdf):
    num_db = get_db_number(vertices_labels)
    g = igraph.Graph()
    g.add_vertices(len(vertices_labels))
    g.vs["label"] = vertices_labels #[s.split(":")[0] for s in vertices_labels]
    colors = [""] * len(vertices_labels)
    for i in range(0, len(vertices_labels)):
        colors[i] = get_color_by_label(vertices_labels[i])
    g.vs["color"] = ["white"] * len(vertices_labels)
    g.vs["shape"] = ["rect" if c == "white" else "circle" for c in g.vs["color"]] # vertices_shapes
#    g.vs["label_dist"] = [0] * len(vertices_labels)
    g.vs["size"] = [30 if sh == "circle" else 30 for sh in g.vs["shape"]]
    g.vs["size2"] = [30 if sh == "circle" else 100 for sh in g.vs["shape"]]
    g.vs["label_color"] = ["black" if c == "white" else c for c in colors] #["white" if sh == 'circle' else "black" for sh in g.vs["shape"]]
    g.vs["label_size"] = [12 if sh == 'circle' else 10 for sh in g.vs["shape"]]
    g.vs["label_degree"] = [45] * len(vertices_labels)
    g.add_edges(edges)
    layout = g.layout_auto() #layout_reingold_tilford(root = get_db_number(vertices_labels)) #(hgap = 1)
    layout.scale(10)
    igraph.plot(g, output_pdf, layout = layout, margin = 75, vertex_frame_width = 0, vertex_label_degree = 20)

def draw_with_igraph_nice(vertices_labels, edges, output_pdf):
    num_db = get_db_number(vertices_labels)
    g = igraph.Graph()
    g.add_vertices(len(vertices_labels))
    g.vs["label"] = [s.split(":")[0] for s in vertices_labels]
    colors = [""] * len(vertices_labels)
    for i in range(0, len(vertices_labels)):
        colors[i] = get_color_by_label(vertices_labels[i])
    g.vs["color"] = colors
    g.vs["shape"] = ["rect" if c == "white" else "circle" for c in g.vs["color"]] # vertices_shapes
#    g.vs["label_dist"] = [0] * len(vertices_labels)
    g.vs["size"] = [30 if sh == "circle" else 30 for sh in g.vs["shape"]]
    g.vs["size2"] = [30 if sh == "circle" else 100 for sh in g.vs["shape"]]
    g.vs["label_color"] = ["white" if sh == 'circle' else "black" for sh in g.vs["shape"]]
    g.vs["label_size"] = [12 if sh == 'circle' else 10 for sh in g.vs["shape"]]
    g.vs["label_degree"] = [45] * len(vertices_labels)
    g.add_edges(edges)
    layout = g.layout_fruchterman_reingold() #(hgap = 1)
    layout.scale(10)
    igraph.plot(g, output_pdf, layout = layout, margin = 75, vertex_frame_width = 0)

def draw_with_graph_tools(vertices_labels, edges, output_pdf):
    g = Graph()
    for i in range(0, len(vertices_labels)):
        g.add_vertex()
    for e in edges:
        g.add_edge(e[0], e[1])
    edge_weight = g.new_edge_property("int")
    for e in g.edges():
        edge_weight[e] = 100
    pos = fruchterman_reingold_layout(g, r = 5, a = 4, circular = True) #(g, eweight = edge_weight, C = 1, p = 7)
    #pos = arf_layout(g, d = 10, max_iter=0)
    vertex_color = g.new_vertex_property("string")
    vertex_label = g.new_vertex_property("string")
    vertex_shape = g.new_vertex_property("string")
    text_position = g.new_vertex_property("int")
    for i in range(0, len(vertices_labels)):
        vertex_color[g.vertex(i)] = get_color_by_label(vertices_labels[i])
        vertex_label[g.vertex(i)] = vertices_labels[i]
        vertex_shape[g.vertex(i)] = get_shape_by_label(vertices_labels[i])
        text_position[g.vertex(i)] = get_text_position(vertices_labels[i])
    g.vertex_properties["vertex_color"] = vertex_color
    g.vertex_properties["vertex_text"] = vertex_label
    #graphviz_draw(g, pos, layout = "arf", overlap = 'scale', output = output_pdf, vcolor = vertex_color,
    #              vprops = {"vertex_text" : vertex_label})
    graph_draw(g, pos, vertex_text = vertex_label, vertex_shape = vertex_shape,
               vertex_fill_color = vertex_color, vertex_font_size = 3, vertex_text_position = text_position,
               vertex_pen_width = 0,
               output_size=(200, 200), output=output_pdf)

def main():
    graph_file = sys.argv[1]
    if not os.path.exists(graph_file):
        print "ERROR: Graph file " + graph_file + " was not found"
        sys.exit(1)
    ifhandler = open(graph_file, "r")
    lines = ifhandler.readlines()
    edge_index = 0
    vertices_colors = []
    vertices_labels = []
    vertices_shapes = []
    for i in range(1, len(lines)):
        if lines[i].strip() == 'Edges:':
            edge_index = i + 1
            break
        splits = lines[i].strip().split()
        vertices_labels.append(splits[1])
        vertices_colors.append(splits[2])
        vertices_shapes.append(splits[3])
    print str(len(vertices_labels)) + " vertices were extracted"
    edges = []
    for i in range(edge_index, len(lines)):
        splits = lines[i].strip().split()
        edges.append((int(splits[0]), int(splits[1])))
    print str(len(edges)) + " edges were extracted"

    output_dir = "" #"graphs"
    output_pdf = os.path.join(output_dir, ntpath.basename(graph_file).split(".")[0] + ".pdf")
    draw_with_igraph_nice(vertices_labels, edges, output_pdf)
    #os.system("open " + output_pdf)
    #draw_with_graph_tools(vertices_labels, edges, output_pdf)

if __name__ == "__main__":
    main()