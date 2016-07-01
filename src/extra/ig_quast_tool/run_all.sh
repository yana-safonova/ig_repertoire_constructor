#!/bin/bash


for i in 1 2 3 4 5 6 7 8 9
do
    # bash make_reference.sh /Nancy/data/input/ImmunoSeq/ibh_datasets/age_datasets/raw_reads/merged_reads/age_ig_s${i}_R12.fastq AGE${i}
    bash make_reference.sh AGE${i}/input_reads.fa.gz AGE_NEW_${i}
    bash run_quast.sh AGE_NEW_${i}
done
