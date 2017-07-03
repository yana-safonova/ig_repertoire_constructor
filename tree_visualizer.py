#! /usr/bin/env python2

import argparse
import subprocess
from os import listdir, rmdir, mkdir
from shutil import rmtree
from math import log 
import os



def abundance_to_size(n):
	coeff = 1
	if n < 10: 
		coeff = 0.8
	elif n < 100:
		coeff = 1.6
	elif n < 1000:
		coeff = 2.4
	else:
		coeff = 3.2

	coeff = log(n+3, 4)
	return coeff*1.2, 1.2 #coeff*0.8

def draw_tree(antevolo_res_dir, tree_name, output_dir):
	tree_file = os.path.join(antevolo_res_dir, 'clonal_trees/', tree_name)
	print "drawing "+tree_file
	vertices_file = os.path.join(antevolo_res_dir,'clonal_trees_vertices/', tree_name)
	vertex_to_depths = {}
	depth_to_vertices = {}
	edges = []
	max_depth = 0
	vertices = {}
	clones = {}
	passed = set()
	with open(vertices_file) as inp:
		inp.readline()
		for st in inp:
			arr = st.split()
			clone_num = int(arr[0])
			clone_name = arr[1]
			clone_AA_seq = arr[3]
			clone_productive = int(arr[2])
			clone_left_anchor_AA = arr[5]
			clone_right_anchor_AA = arr[6]
			clones[clone_num] = [clone_productive, clone_AA_seq, clone_left_anchor_AA, clone_right_anchor_AA]
			clone_abundance = int(arr[7])

			if clone_productive:
				clone_shape = 'ellipse'
			else:
				clone_shape = 'box'

			#clone_abundance = int(clone_name.split('_')[-1].split('|')[0])
			clone_width, clone_height = abundance_to_size(clone_abundance)
			if clone_name.split('_')[0] == 'fake':
				clone_color = 'magenta'
			else:
				clone_color = 'cyan'
			vertices[clone_num] = ''.join(['[label=',"\"" + str(clone_num)+'_'+clone_left_anchor_AA+clone_right_anchor_AA + "\"", 
									 	   ', fixedsize=true, style=filled, fillcolor=', clone_color, 
										   ', shape=', clone_shape, 
										   ' width=', str(clone_width), ' height=', str(clone_height), ']'])

	with open(tree_file) as inp:
		inp.readline()
		for st in inp:
			arr = st.split()
			src_num = int(arr[0])
			dst_num = int(arr[1])
			edge_type = arr[6]
			length = int(arr[8])
			src_depth = int(arr[4])
			dst_depth = int(arr[5])

			vertex_to_depths[src_num] = src_depth
			vertex_to_depths[dst_num] = dst_depth
			
			depth_to_vertices.setdefault(src_depth, set())
			depth_to_vertices[src_depth].add(src_num)
			depth_to_vertices.setdefault(dst_depth, set())
			depth_to_vertices[dst_depth].add(dst_num)
			
			edges.append([src_num, dst_num, edge_type, src_depth, dst_depth])	
			max_depth = max(max_depth, dst_depth)

			passed.add(src_num)
			passed.add(dst_num)

	fake_vertices = ['Depth_'+str(i) for i in xrange(max_depth+1)]
	DOT_OUTPUT_FILE_NAME = os.path.join(output_dir, tree_file.split('/')[-1]+".dot")
	with open(DOT_OUTPUT_FILE_NAME, 'w') as otp:
		otp.write("digraph "+'tree'+' {\n')
		otp.write("\tranksep=equally;")
		otp.write(''.join( ["\t{\n\t\tnode [shape=box]\n\t\t\n\t\t", ' -> '.join(fake_vertices), ";\n\t}\n\n"] ))
		for depth in depth_to_vertices:
			otp.write(''.join( ["\t{ rank = same;\n \t\t", 
								"Depth_"+str(depth)+"; ", 
								"; ".join(["\""+str(num)+"\"" for num in depth_to_vertices[depth]]),
								";\n\t};\n"] ))
		
		for v in vertices:
			if v in passed:
				otp.write(''.join(["\t","\""+str(v)+"\"", vertices[v],";\n"]))

		for edge in edges:
			src_num, dst_num, edge_type, src_depth, dst_depth = edge
			if edge_type == 'undirected' and src_depth == 0:
				continue

			if clones[src_num][1] == clones[dst_num][1]:
				edge_style = 'dotted'
			else:
				edge_style = 'filled'

			if clones[src_num][2] != clones[dst_num][2] or clones[src_num][3] != clones[dst_num][3]:
				edge_color = 'red'
			else:
				edge_color = 'black'

			dir_attr = ''
			if edge_type == 'reverse_directed':
				dir_attr = ', dir=both, arrowhead=none'
				src_num, dst_num = dst_num, src_num
			otp.write(''.join(["\t","\""+str(src_num)+"\"", " -> ", "\""+str(dst_num)+"\"",
				 " [color=", edge_color, ", style=", edge_style, dir_attr, "];\n"]))

		otp.write("}\n")
	subprocess.call(['dot', '-Tpdf', '-O', DOT_OUTPUT_FILE_NAME])


def main():
	parser = argparse.ArgumentParser()

	parser.add_argument("-i", "--input", dest = 'input', help="input dir with AntEvolo results", required=True)
	parser.add_argument("-o", "--output", dest = 'output', help="output dir", required=True)
	parser.add_argument("-s", "--strategy", dest = 'strategy',\
	 help="'single' for specific tree (then -n TREE_FILE_NAME, 'topk' for a number of top-sized trees (then -k NUMBER_OF_TREES)", choices=['single', 'topk'] ,required=True)
	parser.add_argument("-k", "--trees_num", dest = 'k', type=int, help="number of top-size trees to draw")
	parser.add_argument("-n", "--name", dest = 'name', help="file name")
	args = parser.parse_args()

	if args.strategy == 'topk':
		trees = listdir(os.path.join(args.input, "clonal_trees/"))
		trees.sort(key=lambda x: int(x.split('.')[-2].split('_')[-1]))
		#for tree_file in trees[-int(args.k):]:
		#	print tree_file.split('_')[-1].split('.')[-2]
		
		try:
			listdir(args.output)
		except OSError:
			mkdir(args.output)
		for tree_name in trees[-int(args.k):]:
			draw_tree(args.input, tree_name, args.output)
	else:
		draw_tree(args.input, args.name, args.output)

if __name__ == "__main__":
	main()