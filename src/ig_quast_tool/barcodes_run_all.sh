#!/bin/bash

DIR="/johnny/data/input/Ig/ibh_datasets/age_datasets/raw_reads/merged_reads/gzipped"
OUT_DIR = "/home/ashlemov/barcodes_agressive/"

mkdir ${OUT_DIR}

for i in `seq 1 9`; do
    ./barcode_cleaner_nojoin.py -m 5 -b ${OUT_DIR}/age_${i}_bad.fq.gz \
        --rcm $barcodes_${i}.rcm \
        $DIR/age_ig_s${i}_R12_new.fastq.gz \
        ${OUT_DIR}/age_${i}_good.fq.gz &
done

join
