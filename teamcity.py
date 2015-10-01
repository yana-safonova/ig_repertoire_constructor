#!/usr/bin/env python

import os.path
import sys
import shutil
import json

class Params:
    reads_filename = '/Johnny/data/input/Ig/ibh_datasets/age_datasets/raw_reads/cropped_reads/age_ig_s3_R12_raw.cropped.fastq'
    assembled_barcodes_filename = '/Johnny/data/input/Ig/ibh_datasets/age_datasets/assembled_barcodes/cropped_barcodes/age_ig_s3_R12_assembled.cropped.fastq'
    assembled_barcodes_rcm = '/Johnny/data/input/Ig/ibh_datasets/age_datasets/assembled_barcodes/barcodes/rcms_tau_3/barcodes_3.rcm'
    assembly_dir = 'igrepcon_output'
    metrics_dir = 'cmp_output'

def run_barcode_quast():
    cmd = './src/barcode_metrics/barcode_quast.py --Bc ' + Params.assembled_barcodes_filename + \
          ' --Br ' + Params.assembled_barcodes_rcm + \
          ' -c ' + os.path.join(Params.assembly_dir, 'final_repertoire_large.fa') + \
          ' -r ' + os.path.join(Params.assembly_dir, 'final_repertoire.rcm') + \
          ' --threads-num 8 --out ' + Params.metrics_dir
    exit_code = os.system(cmd)
    if exit_code != 0:
        print('BarcodeQuast finished abnormally with exit code ' + str(exit_code))
        sys.exit(5)

def assess_barcode_quast_results():
    metrics_file = os.path.join(Params.metrics_dir, 'metrics.json')
    handler = open(metrics_file)
    metrics_data = json.load(handler)
    handler.close()
    good_barcodes = metrics_data['general_metrics']['#good barcodes']
    bad_barcodes = metrics_data['general_metrics']['#bad barcodes']
    good_barcodes_rate = good_barcodes * 100.0 / (good_barcodes + bad_barcodes)
    if good_barcodes_rate < 85.0:
        print('Good barcodes rate is too low, ' + str(good_barcodes_rate))
        sys.exit(6)
    barcodes_on_dist_0 = metrics_data['general_metrics']['Barcodes on distance 0']
    barcodes_number = metrics_data['individual_metrics']['barcodes']['#clusters']
    data_clusters_number = metrics_data['individual_metrics']['data']['#clusters']
    barcodes_representation = barcodes_on_dist_0 * 100.0 / barcodes_number
    if barcodes_representation < 70.0:
        print('Barcodes representation is too low, ' + str(barcodes_representation))
        sys.exit(7)

def run_test(compile=False):
    if compile:
        if os.path.exists('build'):
            shutil.rmtree('build', True)
        if os.path.exists(Params.assembly_dir):
            shutil.rmtree(Params.assembly_dir)
        if os.path.exists(Params.metrics_dir):
            shutil.rmtree(Params.metrics_dir)
        exit_code = os.system('./prepare_cfg')
        if exit_code != 0:
            print('Preparing configuration files finished abnormally with exit code ' + str(exit_code))
            sys.exit(1)
        exit_code = os.system('make -j 8')
        if exit_code != 0:
            print('Compilation finished abnormally with exit code ' + str(exit_code))
            sys.exit(2)
    cmd = "./ig_repertoire_constructor.py -s %s -o %s -t 8 -C heavy --debug" % (Params.reads_filename, Params.assembly_dir)
    exit_code = os.system(cmd)
    if exit_code != 0:
        print('IgRepCon finished abnormally with exit code ' + str(exit_code))
        sys.exit(3)
    os.system('chmod -R 777 ' + Params.assembly_dir)
    log = open(os.path.join(Params.assembly_dir, 'ig_repertoire_constructor.log')).readlines()
    if len(log) < 2 or log[-2].strip() != 'Thank you for using IgRepertoireConstructor!':
        print('Something goes wrong in IgRepertoireConstructor run')
        sys.exit(4)
    run_barcode_quast()
    assess_barcode_quast_results()

if __name__ == '__main__':
    run_test()
