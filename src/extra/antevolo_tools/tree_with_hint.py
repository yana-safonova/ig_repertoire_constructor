#! /usr/bin/env python2

import argparse
import subprocess
from os import listdir

parser = argparse.ArgumentParser()

parser.add_argument("-t", "--tree_old", help="input file with old clonal tree")
parser.add_argument("-r", "--reads_old", help="input file with old reads")
parser.add_argument("-T", "--tree_new", help="input file with new clonal tree")
parser.add_argument("-R", "--reads_new", help="input file with new reads")
parser.add_argument("-o", "--output", help="output dir")
args = parser.parse_args()


def compute_name_to_read_map(reads_file):
	with open(reads_file) as inp_reads:
		name_to_read = {}
		name = ''
		for st in inp_reads:
			if st[0] == '>':
				if name in name_to_read and len(name_to_read[name]) != 0:
					name_to_read[name] = ''.join(name_to_read[name])
				else:
					name_to_read.pop(name, 1)
				name = st.strip()[1:]
				name_to_read[name] = []
			else:
				name_to_read[name].append(st.strip())
		if name in name_to_read and len(name_to_read[name]) != 0:
			name_to_read[name] = ''.join(name_to_read[name])
		else:
			name_to_read.pop(name, 1)
	return name_to_read

def compute_tree(tree_file, name_to_read):
	tree = {}
	with open(tree_file) as inp:
		inp.readline()
		for st in inp:
			arr = st.split()
			tree.setdefault(name_to_read[arr[2]], [[],''])
			tree.setdefault(name_to_read[arr[3]], [[],''])
			tree[name_to_read[arr[2]]][0].append(name_to_read[arr[3]])
			tree[name_to_read[arr[3]]][1] = name_to_read[arr[2]]
	return tree
		
def compute_reads_anc_relations(tree, read, root, name_to_read, read_to_read_anc_relations):
	for child in tree[read][0]:
		if tree[child][1] == root:
			read_to_read_anc_relations[root][child] = 'parent'
		else:
			read_to_read_anc_relations[root][child] = 'ancestor'
		compute_reads_anc_relations(tree, child, root, name_to_read, read_to_read_anc_relations)
			
def draw_tree_with_hint(new_tree_file, read_to_read_anc_relations, name_to_read, output_dir, old_tree_file):
	print "drawing "+new_tree_file
	vertex_to_depths = {}
	depth_to_vertices = {}
	edges = []
	max_depth = 0
	vertices = {}
	with open(new_tree_file) as inp:
		inp.readline()
		for st in inp:
			arr = st.split()
			src_name = arr[2]
			dst_name = arr[3]
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
			
			#if src_num == 52:
			#	print read_to_read_anc_relations[name_to_read[src_name]], '\n\n\n'
			

			if name_to_read[src_name] not in read_to_read_anc_relations \
			  or name_to_read[dst_name] not in read_to_read_anc_relations[name_to_read[src_name]]: 
				if name_to_read[dst_name] not in read_to_read_anc_relations:
					vertices[dst_num] = ''
				else:
					vertices[dst_num] = '[style=filled, color=yellow]'
			elif read_to_read_anc_relations[name_to_read[src_name]][name_to_read[dst_name]] == 'parent':
				vertices[dst_num] = '[style=filled, color=red]'
			elif read_to_read_anc_relations[name_to_read[src_name]][name_to_read[dst_name]] == 'ancestor':
				vertices[dst_num] = '[style=filled, color=orange]'

			if name_to_read[src_name] in read_to_read_anc_relations and src_num not in vertices:
				vertices[src_num] = '[style=filled, color=violet]'




	fake_vertices = ['Depth_'+str(i) for i in xrange(max_depth+1)]

	DOT_OUTPUT_FILE_NAME = output_dir+"/"+new_tree_file.split('/')[-1] + "__from__" + old_tree_file.split('/')[-1] + ".dot"
	with open(DOT_OUTPUT_FILE_NAME, 'w') as otp:
		otp.write("digraph "+'tree'+' {\n')
		otp.write(''.join( ["\t{\n\t\tnode [shape=box]\n\t\t\n\t\t", ' -> '.join(fake_vertices), ";\n\t}\n\n"] ))
		
		
	   
		
		
		for depth in depth_to_vertices:
			otp.write(''.join( ["\t{ rank = same;\n \t\t", 
								"Depth_"+str(depth)+"; ", 
								"; ".join(["\""+str(num)+"\"" for num in depth_to_vertices[depth]]),
								";\n\t};\n"] ))
			
		for v in vertices:
			otp.write(''.join(["\t","\""+str(v)+"\"", vertices[v], ";\n"]))
		
		
		for edge in edges:
			src_num, dst_num, edge_type, src_depth, dst_depth = edge
			if edge_type == 'directed':
				otp.write(''.join(["\t","\""+str(src_num)+"\"", " -> ", "\""+str(dst_num)+"\"", " [color=black];\n"]))
			elif edge_type == "undirected":
				otp.write(''.join(["\t","\""+str(src_num)+"\"", " -> ", "\""+str(dst_num)+"\"", " [color=blue];\n"]))
		otp.write("}\n")

	draw = False
	for v in vertices:
		if vertices[v] != '':
			draw = True
	if draw:
		subprocess.call(['dot', '-Tpdf', '-O', DOT_OUTPUT_FILE_NAME])
	else:
		print "\tno clones presented in both datasets"
		subprocess.call(['touch', DOT_OUTPUT_FILE_NAME+'.empty.pdf'])

def tree_comparing(old_reads_file, old_tree_file, new_reads_file, new_tree_file, output_dir):
	old_name_to_read = compute_name_to_read_map(old_reads_file)
	new_name_to_read = compute_name_to_read_map(new_reads_file)
	old_tree = compute_tree(old_tree_file, old_name_to_read)
	read_to_read_anc_relations = {}
	for root in old_tree:
		read_to_read_anc_relations[root] = {}
		compute_reads_anc_relations(old_tree, root, root, old_name_to_read, read_to_read_anc_relations)
	#print old_tree[ROOT_SEQ]
	draw_tree_with_hint(new_tree_file, read_to_read_anc_relations, new_name_to_read, output_dir, old_tree_file)


def compare_reads_sets(old_reads_file, new_reads_file, old_tree_file, new_tree_file):
	old_name_to_read = compute_name_to_read_map(old_reads_file)
	new_name_to_read = compute_name_to_read_map(new_reads_file)

	old_reads_set = set()
	for read in old_name_to_read.values():
		old_reads_set.add(read)
	new_reads_set = set()
	for read in new_name_to_read.values():
		new_reads_set.add(read)
	i = 0
	j = 0
	k = 0
	for read in new_reads_set:
		if read in old_reads_set:
			i += 1

	new_tree = compute_tree(new_tree_file, new_name_to_read)
	old_tree = compute_tree(old_tree_file, old_name_to_read)
	for read in new_tree:
		if read in old_reads_set:
			j += 1
	for read in new_tree:
		if read in old_tree:
			k += 1
	return i, j, k
	


def main():
	tree_comparing(args.reads_old, args.tree_old, args.reads_new, args.tree_new, args.output)

if __name__ == "__main__":
	main()
	'''
	i, j, k = compare_reads_sets(args.reads_old, args.reads_new, args.tree_old, args.tree_new)
	print 'number of reads from 2nd dataset presented in 1st ', i,\
		  'number of reads from 2nd tree presented in 1st dataset', j,\
		  'number of reads from 2nd tree presented in 1st tree', k
	'''