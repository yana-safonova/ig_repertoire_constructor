import os
import sys

convert_command_line = "python convert_to_graph.py"
metis_command_line = "metis-5.1.0/build/Linux-x86_64/programs/ndmetis"
drawer_command_line = "python draw_matrix.py"

def RunGraphConvertion(graph_file):
    err_code = os.system(convert_command_line + " " + graph_file)
    if err_code != 0:
        print "Error: Convertion to GRAPH format was failed"
        return ""
    output_graph = graph_file + ".graph"
    if not os.path.exists(output_graph):
        print "Error: Output graph was not created"
        return ""
    print "Metis graph was written to " + output_graph
    return output_graph

def RunMetis(graph_file):
    err_code = os.system(metis_command_line + " " + graph_file)
    if err_code != 0:
        print "Error: ND Metis finished abnormally"
        return ""
    output_perm = graph_file + ".iperm"
    if not os.path.exists(output_perm):
        print "Error: Permutation was not created"
        return ""
    print "Metis permutation was written to " + output_perm
    return output_perm

def RunDrawer(graph, perm, mode, output_dir):
    err_code = os.system(drawer_command_line + " " + graph + " " + perm + " " + mode)
    if err_code != 0:
        print "Error: Drawer finished abnormally"
        return ""
    output_ppm = os.path.join(output_dir, graph + "_" + mode + ".ppm")
    if not os.path.exists(output_ppm):
        print "Error: PPM file was not created"
        return ""
    print "Output figure was written to " + output_ppm
    return output_ppm

input_dir = os.path.abspath(sys.argv[1])
output_dir = os.path.abspath(sys.argv[2])
input_files = os.listdir(input_dir)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for f in input_files:
    graph_file = os.path.join(input_dir, f)
    print "==========================================="
    print graph_file 
    output_graph = RunGraphConvertion(graph_file)    
    if graph_file != "":
        output_perm = RunMetis(output_graph)
        if output_perm != "":
            RunDrawer(graph_file, output_perm, "weight", output_dir)
            RunDrawer(graph_file, output_perm, "shm", output_dir)
