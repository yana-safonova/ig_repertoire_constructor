#! /usr/bin/env python2

import argparse
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="input .tree file")
parser.add_argument("-o", "--output", help="output file")
args = parser.parse_args()

def main():
	tree = {}
	max_depth = 0
	tree_name = args.input.split('/')[-1].split('.')[-2].split('-')[-1]
	with open(args.input) as inp:
		with open(args.output, 'w') as otp:
			otp.write("digraph "+'tree'+' {\n')
			inp.readline()
			for st in inp:
				arr = st.split()
				src_num = arr[0]
				dst_num = arr[1]
				edge_type = arr[6]
				length = arr[8]
				
				if edge_type == 'directed':
					otp.write(''.join(["\t",src_num, " -> ", dst_num, " [color=black];\n"]))
				elif edge_type == "undirected":
					otp.write(''.join(["\t",src_num, " -> ", dst_num, " [color=blue];\n"]))
			otp.write("}\n")
	subprocess.call(['dot', '-Tpdf', '-O', args.output])

if __name__ == "__main__":
	main()