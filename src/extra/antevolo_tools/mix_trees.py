#! /usr/bin/env python2

import argparse
import subprocess
from os import listdir
from sys import argv
from math import log

parser = argparse.ArgumentParser()

parser.add_argument("-t", "--tree_old", help="input file with old clonal tree")
parser.add_argument("-r", "--reads_old", help="input file with old reads")
parser.add_argument("-T", "--tree_new", help="input file with new clonal tree")
parser.add_argument("-R", "--reads_new", help="input file with new reads")
parser.add_argument("-o", "--output", help="output dir")
args = parser.parse_args()
IGREC_DIR = '/'.join(argv[0].split('/')[:-1])+'/../../../'
THREADS_NUM = '16'

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
			
def draw_colored_tree(mixed_tree_file, output_dir, old_tree_file, new_tree_file):
	print "drawing "+mixed_tree_file
	mixed_tree_name = mixed_tree_file.split('/')[-1]
	vertices_file = '/'.join(mixed_tree_file.split('/')[:-1]) + '/vertices_' + mixed_tree_name
	print vertices_file

	aa_sequences = {}
	with open(vertices_file) as inp:
		inp.readline()
		for st in inp:
			arr = st.split()
			num = int(arr[0])
			aa_sequences[num] = arr[3]


	vertex_to_depths = {}
	depth_to_vertices = {}
	edges = []
	max_depth = 0
	vertices = {}
	with open(mixed_tree_file) as inp:
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
			src_productive = int(arr[11])
			dst_productive = int(arr[12])
			synonymous = int(arr[13])

			vertex_to_depths[src_num] = src_depth
			vertex_to_depths[dst_num] = dst_depth
			
			depth_to_vertices.setdefault(src_depth, set())
			depth_to_vertices[src_depth].add(src_num)
			depth_to_vertices.setdefault(dst_depth, set())
			depth_to_vertices[dst_depth].add(dst_num)
			
			edges.append([src_num, dst_num, edge_type, src_depth, dst_depth, synonymous])	
			max_depth = max(max_depth, dst_depth)
			
			#if src_num == 52:
			#	print read_to_read_anc_relations[name_to_read[src_name]], '\n\n\n'
			
			src_presented_in = src_name.split('|')[-1]
			if src_presented_in == 'old':
				src_color = 'red'
			elif src_presented_in == 'new':
				src_color = 'blue'
			elif src_presented_in == 'both':
				src_color = 'orange'
			dst_presented_in = dst_name.split('|')[-1]
			if dst_presented_in == 'old':
				dst_color = 'red'
			elif dst_presented_in == 'new':
				dst_color = 'blue'
			elif dst_presented_in == 'both':
				dst_color = 'orange'


			if src_productive:
				src_shape = 'ellipse'
			else:
				src_shape = 'box'
			if dst_productive:
				dst_shape = 'ellipse'
			else:
				dst_shape = 'box'

			src_abundance = int(src_name.split('_')[-1].split('|')[0])
			src_width, src_height = abundance_to_size(src_abundance)
			dst_abundance = int(dst_name.split('_')[-1].split('|')[0])
			dst_width, dst_height = abundance_to_size(dst_abundance)

			vertices[src_num] = ''.join(['[fixedsize=true, style=filled, fillcolor=', src_color, ", shape=", src_shape, 
										 " width=", str(src_width), " height=", str(src_height), ']'])
			vertices[dst_num] = ''.join(['[fixedsize=true, style=filled, fillcolor=', dst_color, ", shape=", dst_shape,
										 ' width=', str(dst_width), " height=", str(dst_height), ']'])

			




	fake_vertices = ['Depth_'+str(i) for i in xrange(max_depth+1)]

	DOT_OUTPUT_FILE_NAME = output_dir+"/mixed_"+new_tree_file.split('/')[-1] + "__" + old_tree_file.split('/')[-1] + ".dot"
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
			src_num, dst_num, edge_type, src_depth, dst_depth, synonymous = edge
			if edge_type == "undirected" and src_depth == 0:
				continue

			synonymous = aa_sequences[src_num] == aa_sequences[dst_num]
			if synonymous:
				otp.write(''.join(["\t","\""+str(src_num)+"\"", " -> ", "\""+str(dst_num)+"\"", " [color=magenta];\n"]))
			else:
				otp.write(''.join(["\t","\""+str(src_num)+"\"", " -> ", "\""+str(dst_num)+"\"", " [color=black];\n"]))
		otp.write("}\n")

		'''
		for edge in edges:
			src_num, dst_num, edge_type, src_depth, dst_depth = edge
			if edge_type == 'directed':
				otp.write(''.join(["\t","\""+str(src_num)+"\"", " -> ", "\""+str(dst_num)+"\"", " [color=black];\n"]))
			elif edge_type == "undirected":
				otp.write(''.join(["\t","\""+str(src_num)+"\"", " -> ", "\""+str(dst_num)+"\"", " [color=blue];\n"]))
		otp.write("}\n")
		'''
	draw = False
	for v in vertices:
		if vertices[v] != '':
			draw = True
	if draw:
		subprocess.call(['dot', '-Tpdf', '-O', DOT_OUTPUT_FILE_NAME])
	else:
		print "\tno clones presented in both datasets"
		subprocess.call(['touch', DOT_OUTPUT_FILE_NAME+'.empty.pdf'])

def mix_trees(old_reads_file, old_tree_file, new_reads_file, new_tree_file):
	
	antevolo_res_dir = '/temp_antevolo_res_'+new_tree_file.split('/')[-1] + "__" + old_tree_file.split('/')[-1] + '/'

	old_name_to_read = compute_name_to_read_map(old_reads_file)
	new_name_to_read = compute_name_to_read_map(new_reads_file)


	reads = {}
	old_names = set()
	new_names = set()
	with open(old_tree_file) as inp:
		inp.readline()
		for st in inp:
			arr = st.split()
			old_names.add(arr[2])
			old_names.add(arr[3])
	with open(new_tree_file) as inp:
		inp.readline()
		for st in inp:
			arr = st.split()
			new_names.add(arr[2])
			new_names.add(arr[3])


	for name in old_names:
		read = old_name_to_read[name]
		reads.setdefault(read, {'old': '', 'new' : ''})
		reads[read]['old'] = name
	for name in new_names:
		read = new_name_to_read[name]
		reads.setdefault(read, {'old': '', 'new' : ''})
		reads[read]['new'] = name

	# prepare the antevolo input
	try :
		subprocess.check_output(['mkdir', args.output+antevolo_res_dir])
	except Exception:
		subprocess.check_output(['rm', '-rf',  args.output+antevolo_res_dir])

	both = 0
	old = 0
	new = 0
	with open(args.output+'/mixed_reads.fa', 'w') as otp_reads:
		for read in reads:
			if reads[read]['old'] == '' and reads[read]['new'] == '':
				raise Exception("Error: read is not presented in the datasets")
			if reads[read]['old'] != '' and reads[read]['new'] == '':
				name = ''.join(['>', reads[read]['old'], '|old', '\n'])
				otp_reads.write(name)
				otp_reads.write(read+'\n')
				old += 1
			if reads[read]['old'] == '' and reads[read]['new'] != '':
				name = ''.join(['>', reads[read]['new'], '|new', '\n'])
				otp_reads.write(name)
				otp_reads.write(read+'\n')
				new += 1
			if reads[read]['old'] != '' and reads[read]['new'] != '':
				name = ''.join(['>', reads[read]['old'], '|both', '\n'])
				otp_reads.write(name)
				otp_reads.write(read+'\n')
				both += 1

	print '\n\n\n', 'both: ', both, 'new: ', new, 'old: ', old, '\n\n\n'
	subprocess.call([IGREC_DIR+'antevolo.py', 
					 '-i', args.output+'/mixed_reads.fa',
					 '-o', args.output+antevolo_res_dir,
					 '-t', THREADS_NUM,
					 '-l', 'IG'])
	
	mixed_trees = listdir(args.output+antevolo_res_dir+'clonal_trees/')
	mixed_trees.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
	subprocess.call(['mv', args.output+antevolo_res_dir+'/clonal_trees/'+mixed_trees[-1], args.output+'/'])
	subprocess.call(['mv', args.output+antevolo_res_dir+'/clonal_trees_vertices/'+mixed_trees[-1], args.output+'/vertices_'+mixed_trees[-1]])
	subprocess.call(['rm', '-rf', args.output+antevolo_res_dir])
	subprocess.call(['rm', args.output+'mixed_reads.fa'])

	draw_colored_tree(args.output+'/'+mixed_trees[-1], args.output+'/', old_tree_file, new_tree_file)
	'''
	old_tree = compute_tree(old_tree_file, old_name_to_read)
	read_to_read_anc_relations = {}
	for root in old_tree:
		read_to_read_anc_relations[root] = {}
		compute_reads_anc_relations(old_tree, root, root, old_name_to_read, read_to_read_anc_relations)
	#print old_tree[ROOT_SEQ]
	draw_tree_with_hint(new_tree_file, read_to_read_anc_relations, new_name_to_read, output_dir, old_tree_file)
	'''

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
	mix_trees(args.reads_old, args.tree_old, args.reads_new, args.tree_new)

if __name__ == "__main__":
	main()
	'''
	i, j, k = compare_reads_sets(args.reads_old, args.reads_new, args.tree_old, args.tree_new)
	print 'number of reads from 2nd dataset presented in 1st ', i,\
		  'number of reads from 2nd tree presented in 1st dataset', j,\
		  'number of reads from 2nd tree presented in 1st tree', k
	'''