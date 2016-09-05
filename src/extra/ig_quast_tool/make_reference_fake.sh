#!/bin/bash


INPUT=$1
OUTPUT=$2

IGREC_DIR=/home/ashlemov/Git/ig_repertoire_constructor

BCLEANER=${IGREC_DIR}/src/extra/ig_quast_tool/barcode_cleaner_nojoin.py
CFINDER="${IGREC_DIR}/build/release/bin/ig_consensus_finder -H"
CIC="${IGREC_DIR}/py/ig_compress_equal_clusters.py "

${BCLEANER} ${OUTPUT}/cleaned_reads.fa.gz ${OUTPUT}/input_fake.fa.gz -r ${OUTPUT}/repertoire_fake.rcm --tau=0 -d100500 --lengths ${OUTPUT}/lengths.txt
