#!/bin/bash

DIR="/johnny/data/input/Ig/ibh_datasets/age_datasets/raw_reads/merged_reads/gzipped"
OUT_DIR="/home/ashlemov/barcodes_agressive/"

mkdir ${OUT_DIR}

for i in `seq 1 9`; do
    ../../build/release/bin/ig_consensus_finder -H \
        -i ${OUT_DIR}/age_${i}_good.fq.gz \
        -R ${OUT_DIR}/barcodes_${i}.rcm \
        -o ${OUT_DIR}/repertoire_${i}.fa \
        -t 2 &
done

wait
