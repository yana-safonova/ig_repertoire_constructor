from Bio import SeqIO
from repertoire_structs import Repertoire, Cluster
from collections import defaultdict
import sys
import os.path

def read_repertoire(clusters_file, rcm_file, name):
    clusters = None
    if clusters_file.endswith('.fa') or clusters_file.endswith('.fasta'):
        clusters = read_clusters_fa(clusters_file)
    else:
        clusters = read_migec(clusters_file)
    cluster_reads, read_clusters = None, None
    if rcm_file:
        cluster_reads, read_clusters = read_rcm(rcm_file)
    return Repertoire(clusters_file, rcm_file, clusters, cluster_reads, read_clusters, name)

def read_repertoires(barcode_clusters, barcode_rcm, data_clusters, data_rcm):
    barcode_repertoire = read_repertoire(barcode_clusters, barcode_rcm, 'assembled_barcodes')
    name = os.path.basename(data_clusters).split('.', 1)[0]
    data_repertoire = read_repertoire(data_clusters, data_rcm, name)
    return barcode_repertoire, data_repertoire

def cut_and_get_abundance(id_string):
    abundance = 1
    if 'abundance' in id_string:
        abundance = int(id_string.split('_')[-1].split(':')[1])
        id_string = '_'.join(id_string.split('_')[:-1])
    return abundance, id_string

def parse_cluster_fa_id(id_string):
    abundance, id_string = cut_and_get_abundance(id_string)
    fields = id_string.split('___')
    if len(fields) < 4:
        print "ERROR: wrong format of cluster id " + string_ids + ", must be in format cluster___id___size___sz"
        sys.exit(-1)
    return int(fields[1]), int(fields[3]), abundance

def parse_migec_id(id_string):
    abundance, id_string = cut_and_get_abundance(id_string)
    fields = id_string.split(':')[-1].split('_')
    if len(fields) != 2:
        print "ERROR: wrong format of MiGEC ids, must be in format @MIG_UMI:BARCODE:size_id"
        sys.exit(-1)
    return int(fields[1]), int(fields[0]), abundance

def parse_seq_id(id_string):
    if '___' in id_string:
        cluster_id, cluster_size, cluster_abundance = parse_cluster_fa_id(id_string)
    else:
        cluster_id, cluster_size, cluster_abundance = parse_migec_id(id_string)
    return cluster_id, cluster_size, cluster_abundance

def read_clusters_fa(fasta):
    clusters = {}
    fasta_records = SeqIO.parse(open(fasta), 'fasta')
    for rec in fasta_records:
        cluster_id, cluster_size, cluster_abundance = parse_seq_id(rec.id)
        clusters[cluster_id] = Cluster(cluster_size, rec.seq, cluster_abundance)
    return clusters

def read_migec(migec_filename):
    clusters = {}
    fastq_records = SeqIO.parse(open(migec_filename), 'fastq')
    for rec in fastq_records:
        cluster_id, cluster_size, cluster_abundance = parse_seq_id(rec.id)
        clusters[cluster_id] = Cluster(cluster_size, rec.seq, cluster_abundance)
    return clusters

def read_rcm(rcm):
    cluster_reads = defaultdict(set)
    read_clusters = {}
    for l in open(rcm):
        fields = l.strip().split('\t')
        read_id = fields[0]
        if len(fields) != 2:
            print "ERROR: wrong line in rcm: " + l.strip() + ", must be in format read_id\\tcluster_no"
            sys.exit(1)
        cluster_id = fields[1]
        read_clusters[read_id] = cluster_id
        cluster_reads[cluster_id].add(read_id)
    return cluster_reads, read_clusters

def read_cluster_numbers_to_ids(filtered_clusters_fa):
    cluster_num_to_id = {}
    fasta_records = SeqIO.parse(open(filtered_clusters_fa), 'fasta')
    for i, rec in enumerate(fasta_records):
        cluster_id, cluster_size, abundance = parse_seq_id(rec.id)
        cluster_num_to_id[i] = cluster_id
    return cluster_num_to_id

def read_cluster_matches(matches_filename, repertoire, this_cluster_num_to_ids, other_cluster_num_to_ids, tau):
    cluster_matches = {}
    handler = open(matches_filename)
    for l in open(matches_filename):
        for this_cluster_no, l in enumerate(handler):
            fields = l.strip().split(' ')
            if not fields[0]:
                continue
            score = -int(fields[0])
            if score > tau:
                continue
            cluster_nums = (int(i) for i in fields[1:])
            this_cluster_id = this_cluster_num_to_ids[this_cluster_no]
            cluster_matches[this_cluster_id] = [set(), score]
            for other_cluster_no in cluster_nums:
                other_cluster_id = other_cluster_num_to_ids[other_cluster_no]
                cluster_matches[this_cluster_id][0].add(other_cluster_id)
    return cluster_matches
