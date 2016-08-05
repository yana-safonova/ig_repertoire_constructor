#! /usr/bin/env python2

import argparse
import subprocess
from os import listdir

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="input dir with clonal trees")
parser.add_argument("-o", "--output", help="output dir")
parser.add_argument("-s", "--strategy", help="dot/neato/oth")
parser.add_argument("-k", "--trees_num", help="number of top-size trees to draw")
args = parser.parse_args()

def draw_tree(tree_file):
	print "drawing "+tree_file
	vertex_to_depths = {}
	depth_to_vertices = {}
	edges = []
	max_depth = 0
	vertices = {}
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

			vertices[dst_num] = 'yellow'


	fake_vertices = ['Depth_'+str(i) for i in xrange(max_depth+1)]


	DOT_OUTPUT_FILE_NAME = args.output+"/"+tree_file.split('/')[-1]+".dot"
	with open(DOT_OUTPUT_FILE_NAME, 'w') as otp:
		otp.write("digraph "+'tree'+' {\n')
		otp.write(''.join( ["\t{\n\t\tnode [shape=box]\n\t\t\n\t\t", ' -> '.join(fake_vertices), ";\n\t}\n\n"] ))

		for depth in depth_to_vertices:
			otp.write(''.join( ["\t{ rank = same;\n \t\t", 
								"Depth_"+str(depth)+"; ", 
								"; ".join(["\""+str(num)+"\"" for num in depth_to_vertices[depth]]),
								";\n\t};\n"] ))
		
		#for v in vertices:
		#	otp.write(''.join(["\t","\""+str(v)+"\"", " [style=filled, color=", vertices[v], "];\n"]))

		for edge in edges:
			src_num, dst_num, edge_type, src_depth, dst_depth = edge
			if edge_type == 'undirected' and src_depth == 0:
				continue
			if edge_type == 'directed':
				otp.write(''.join(["\t","\""+str(src_num)+"\"", " -> ", "\""+str(dst_num)+"\"", " [color=black];\n"]))
			elif edge_type == "undirected":
				otp.write(''.join(["\t","\""+str(src_num)+"\"", " -> ", "\""+str(dst_num)+"\"", " [color=blue];\n"]))


		otp.write("}\n")
	subprocess.call([args.strategy, '-Tpdf', '-O', DOT_OUTPUT_FILE_NAME])


def main():
	trees = listdir(args.input)
	trees.sort(key=lambda x: int(x.split('.')[-2].split('_')[-1]))
	for tree_file in trees[-int(args.trees_num):]:
		draw_tree(args.input+'/'+tree_file)

if __name__ == "__main__":
	main()