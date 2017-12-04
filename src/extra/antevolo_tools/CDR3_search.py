#! /usr/bin/env python2

import argparse
import subprocess
from os import listdir

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="input dir with clonal trees")
parser.add_argument("-o", "--output", help="output dir")
parser.add_argument("-c", "--cdr", help="CDR3 string")
args = parser.parse_args()

CDR3_to_find = args.cdr
tau = 3
#CDR3_to_find = 'GCATTCAGCAGCTGGCAGGGGGAAGTTGACTAC' # 9830, 4th 
#CDR3_to_find = 'GTCACACTGTTTGGGTGTGACGAC', 49832, 5th

def hamming_distance(a, b):
	if len(a) != len(b):
		return max(len(a), len(b))
	res = 0
	for i in xrange(len(a)):
		if a[i] != b[i]:
			res += 1
	return res

def search_for_CDR3(tree_file):
	related_trees = set()
	found = False
	#print "drawing "+tree_file
	vertex_to_depths = {}
	depth_to_vertices = {}
	edges = []
	max_depth = 0
	with open(tree_file) as inp:
		inp.readline()
		for st in inp:
			arr = st.split()
			
			src_id = int(arr[0])
			dst_id = int(arr[1])
			'''
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
			'''
			src_CDR3 = arr[14]
			dst_CDR3 = arr[15]

			if hamming_distance(CDR3_to_find, src_CDR3) <= tau or hamming_distance(CDR3_to_find, dst_CDR3) <= tau:
				found = True
				#print '\n\n'+CDR3_to_find+'\n'+src_CDR3+'\n'+dst_CDR3+'\n\n'
	if found:
		related_trees.add(tree_file.split('/')[-1])

	#print '\n'.join(related_trees), '\n\n\n\n' 
	if len(related_trees) != 0:
		print '\n'.join(related_trees)

def main():
	trees = listdir(args.input)
	trees.sort(key=lambda x: int(x.split('.')[-2].split('_')[-1]))
	for tree_file in trees: #[-int(args.trees_num):]:
		search_for_CDR3(args.input+'/'+tree_file)


if __name__ == "__main__":
	main()