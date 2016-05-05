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
import files_utils

def GetClusterIdMap(local_ids, num_previous_clusters):
    id_map = dict()
    shift = 1
    for cur_id in local_ids:
        id_map[cur_id] = num_previous_clusters + shift
        shift += 1
    return id_map

def ParseRCM(filename):        
    local_ids = set()
    id_readname = dict()
    lines = open(filename).readlines()
    for l in lines:
        splits = l.split("\t")
        cur_id = int(splits[1])
        cur_name = splits[0]
        local_ids.add(cur_id)
            
        if cur_id in id_readname:
            id_readname[cur_id].append(cur_name)
        else:
            id_readname[cur_id] = list()
            id_readname[cur_id].append(cur_name)
    return local_ids, id_readname

def AppendToNewRCM(file_handler, id_map, id_readname):
    for var, value in id_readname.items():
        for read in value:
            file_handler.write(read + "\t" + str(id_map[var]) + "\n")

def ParseClusterFa(filename):
    lines = open(filename).readlines()
    parsed_lines = list()
    if len(lines) % 2 != 0:
        print("File " + filename + "is invalid")
        sys.exit(1)
    for i in range(0, len(lines) / 2):
        header = lines[i * 2].strip('\n')
        header = header[1:len(header)]
        seq = lines[i * 2 + 1].strip('\n')
        splits = header.split("___")
        parsed_lines.append([splits[1], splits[3], seq])
    return parsed_lines

def AppendToNewClusterFa(file_handler, parsed_lines, id_map):
    for l in parsed_lines:
        cur_id = l[0]
        cur_size = l[1]
        cur_seq = l[2]
        file_handler.write(">cluster___" + str(id_map[int(cur_id)]) + "___size___" + str(cur_size) + "\n")
        file_handler.write(cur_seq + "\n")

def CombineClusters(rcm_dir, fa_dir, output_fname):
    general_rcm = open(output_fname + ".rcm", "w")
    general_fa = open(output_fname + ".clusters.fa", "w")

    files = [ f for f in os.listdir(rcm_dir) if os.path.isfile(os.path.join(rcm_dir,f)) ]
    num_clusters = 0
    clusters = dict()
    for f in files:
        local_ids, id_readname = ParseRCM(os.path.join(rcm_dir, f))
        id_map = GetClusterIdMap(local_ids, num_clusters)

        AppendToNewRCM(general_rcm, id_map, id_readname)
        
        fa_fname = os.path.join(fa_dir, files_utils.GetFilenameAndExtension(f)[0] + ".clusters.fa")
        if not os.path.exists(fa_fname):
            print("Invalid RCM and FA directories")
            sys.exit(1)
        parsed_fa_lines = ParseClusterFa(fa_fname)
        AppendToNewClusterFa(general_fa, parsed_fa_lines, id_map)

        num_clusters += len(local_ids)
    return num_clusters

def NumberOfClusters(clusters_fa):
    lines = open(clusters_fa).readlines()
    if len(lines) % 2 != 0:
        print("File " + clusters_fa + "is invalid")
        sys.exit(1)
    return len(lines) / 2
    
def PrintMetricsDictToFile(metrics, fname):
    fhandler = open(fname, "w")
    fhandler.write("#constructed clusters\t\t\t\t" + str(metrics["#constructed clusters"]) + "\n")
    fhandler.write("#constructed singletons\t\t\t\t" + str(metrics["#constructed singletons"]) + "\n")
    fhandler.write("max constructed cluster\t\t\t\t" + str(metrics["max constructed cluster"]) + "\n")
    fhandler.write("avg cluster size\t\t\t\t" + str(metrics["avg cluster size"]) + "\n")

def PrintMetricsDictToCSV(metrics, fname):
    fhandler = open(fname, "w")
    fhandler.write("#constructed clusters\t" + str(metrics["#constructed clusters"]) + "\n")
    fhandler.write("#constructed singletons\t" + str(metrics["#constructed singletons"]) + "\n")
    fhandler.write("max constructed cluster\t" + str(metrics["max constructed cluster"]) + "\n")
    fhandler.write("avg cluster size\t" + str(metrics["avg cluster size"]) + "\n")


def PrintMetricsDict(metrics, log):
    log.info("#constructed clusters\t\t\t\t" + str(metrics["#constructed clusters"]))
    log.info("#constructed singletons\t\t\t\t" + str(metrics["#constructed singletons"]))
    log.info("max constructed cluster\t\t\t\t" + str(metrics["max constructed cluster"]))
    log.info("avg cluster size\t\t\t\t" + str(metrics["avg cluster size"]))
