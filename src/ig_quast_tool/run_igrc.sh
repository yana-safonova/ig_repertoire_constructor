#!/bin/bash

OUT_DIR="/home/ashlemov/barcodes_agressive/"

for i in `seq 1 9`; do
    ../../ig_repertoire_constructor.py -t 16 \
        -s ${OUT_DIR}/age_${i}_good.fq.gz \
        -o ${OUT_DIR}/igrc_out_${i} --debug
done
