import os
import sys

prefix_dir = "abpair_test/demultiplexed_pairing_data"
input_dir = os.path.join(prefix_dir, sys.argv[1])
files_list = os.listdir(input_dir)
for fasta_file in files_list:
    if fasta_file[len(fasta_file) - len("fasta"):] != "fasta":
        continue
    fasta_file = os.path.join(input_dir, fasta_file)
    os.system("./clustalw2 -infile=" + fasta_file)
