#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import getopt
import os
import logging
import shutil

from Bio.Seq import Seq, translate
from os.path import isfile, isdir, join
from os import listdir, curdir

from Bio import SeqIO

def GetFilenameAndExtension(filename):
    splits = filename.split('.')
    fname = ""
    for i in range(0, len(splits) - 2):
        fname = fname + splits[i] + "."
    fname = fname + splits[len(splits) - 2]
    ext = splits[len(splits) - 1]   
    return [fname, ext]

class DataFrame:
    data = dict()
    colnames = list()

def PrintDataFrame(data_frame):
    print(data_frame.colnames)
    for name in data_frame.colnames:
        print(name + ":")
        print(data_frame.data[name])

def ReadData(filename):
    if not os.path.exists(filename):
        print("File " + filename + " is not exists")
        sys.exit(1)
    file_handler = open(filename, 'r')
    
    data_frame = DataFrame()
#    print(data_frame)
    first_line = True
    for line in file_handler.readlines():
        splits = line.strip().split()
#        print("Splits: " + str(splits))
        if first_line:
            for i in range(0, len(splits)):
                new_name = "col" + str(i + 1)
                data_frame.colnames.append(new_name)
                data_frame.data[new_name] = list()
                data_frame.data[new_name].append(splits[i])
            #print(data_frame.colnames)
            first_line = False
        for i in range(0, len(splits)):
            #print(str(i) + " " + str(range(0, len(splits))))
            data_frame.data[data_frame.colnames[i]].append(splits[i])
    #print(data_frame.colnames)

    return data_frame

def StrListToInt(str_list):
    for i in range(0, len(str_list)):
        str_list[i] = int(str_list[i])
    return str_list

def StrListToFloat(str_list):
    for i in range(0, len(str_list)):
        str_list[i] = float(str_list[i])
    return str_list

def remove_extension(s):
    arr = s.split('.')
    if len(arr) < 2:
        return s
    s = ''
    for i in range(len(arr) - 2):
        s += arr[i] + '.'
    s += arr[-2]
    return s

def WriteReadInFastqFile(name, seq, qual, reverse, fastq_fname):
    fastq_file = open(fastq_fname, 'a')
    fastq_file.write("@" + name + "\n")
    seq_to_write = seq
    if reverse:
        seq_to_write = str(Seq(seq).reverse_complement())
    fastq_file.write(seq_to_write + "\n")
    fastq_file.write("+\n")
    qual_to_write = qual
    if reverse:
        qual_to_write = qual[::-1]
    fastq_file.write(qual_to_write + "\n")

def FastqToFasta(infile, outfile):
    SeqIO.convert(infile, "fastq", outfile, "fasta")
    return outfile

def WriteListToFile(data_list, fname):
    fhandler = open(fname, "w")
    for item in data_list:
        fhandler.write(str(item) + "\n")
