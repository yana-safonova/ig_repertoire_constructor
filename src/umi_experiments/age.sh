#!/bin/bash
for i in `seq 1 9`;
do
    ./print_umi_graph_stats.py -i "/Johnny/data/input/Ig/ibh_datasets/age_datasets/raw_reads/merged_reads/age_ig_s${i}_R12.fastq" -o "/Johnny/data/input/Ig/ibh_datasets/age_datasets/raw_reads/merged_reads/stats/age_ig_s${i}_R12.stats" -t "/Johnny/data/input/Ig/ibh_datasets/age_datasets/raw_reads/merged_reads/tmp${i}" &> "age_log${i}" &
done

