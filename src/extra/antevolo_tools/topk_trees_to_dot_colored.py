#! /usr/bin/env python2

import argparse
import subprocess
from os import listdir, rmdir, mkdir
from shutil import rmtree
from math import log 

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="input dir with clonal trees")
parser.add_argument("-o", "--output", help="output dir")
parser.add_argument("-s", "--strategy", help="dot/neato/oth")
parser.add_argument("-k", "--trees_num", help="number of top-size trees to draw")
args = parser.parse_args()


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

def draw_tree(tree_file):
	print "drawing "+tree_file
	tree_name = tree_file.split('/')[-1].split('.')[-2]
	vertices_file = args.input + '/../clonal_trees_vertices/' + tree_name + '.tree'
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

			if clone_productive:
				clone_shape = 'ellipse'
			else:
				clone_shape = 'box'

			clone_abundance = int(clone_name.split('_')[-1].split('|')[0])
			#clone_abundance = 1
			clone_width, clone_height = abundance_to_size(clone_abundance)

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

			#vertices[dst_num] = 'yellow'


	fake_vertices = ['Depth_'+str(i) for i in xrange(max_depth+1)]


	DOT_OUTPUT_FILE_NAME = args.output+"/"+tree_file.split('/')[-1]+".dot"
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


			otp.write(''.join(["\t","\""+str(src_num)+"\"", " -> ", "\""+str(dst_num)+"\"",
				 " [color=", edge_color, ", style=", edge_style, "];\n"]))

		otp.write("}\n")
	subprocess.call([args.strategy, '-Tpdf', '-O', DOT_OUTPUT_FILE_NAME])


def main():
	trees = listdir(args.input)
	trees.sort(key=lambda x: int(x.split('.')[-2].split('_')[-1]))
	for tree_file in trees[-int(args.trees_num):]:
		print tree_file.split('_')[-1].split('.')[-2]
	try:
		rmtree(args.output)
	except OSError:
		try:
			rmdir(args.output)
		except OSError:
			pass
	finally:
		mkdir(args.output)
	for tree_file in trees[-int(args.trees_num):]:
		draw_tree(args.input+'/'+tree_file)

if __name__ == "__main__":
	main()