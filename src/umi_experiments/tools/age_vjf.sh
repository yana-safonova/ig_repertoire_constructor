#!/bin/bash
for i in `seq 3 9`;
do
    ./ig_kplus_vj_finder -i "/Johnny/data/input/Ig/ibh_datasets/age_datasets/raw_reads/merged_reads/age_ig_s${i}_R12.fastq" -o "/Johnny/data/input/Ig/ibh_datasets/age_datasets/raw_reads/merged_reads/vdj/${i}" &> "/home/sbankevich/data/age/vjf${i}.log" &
done

