#!/bin/bash

INPUT_DIR=$1
OUTPUT_DIR=$2

GEN_FREQ_DIR='/home/ashlemov/Yandex.Disk/Documents/lab/vers2/src/andy/gen_clusterisation'


for file in ${INPUT_DIR}/age_ig_s?_R12_assembled.cleaned.fastq; do
  echo ${file}
  filename=$(basename "$file")
  i="${filename:8:1}"
  echo $i
  ${GEN_FREQ_DIR}/gen_freq.py -s ${file} -a "${INPUT_DIR}/igblast_output/age_ig_s${i}_igblast_cleaned.align" -o ${OUTPUT_DIR}/${i} 
done
wait
