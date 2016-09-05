#!/bin/bash

for LOCI in IG
do
    bash make_reference_fake.sh /Nancy/data/input/ImmunoSeq/AbVitro/flu_time_course/FV/single_reads_with_UMIs/FV_21.fastq  FLU_FV_21_${LOCI} $LOCI &
    bash make_reference_fake.sh /Nancy/data/input/ImmunoSeq/AbVitro/flu_time_course/FV/single_reads_with_UMIs/FV_22.fastq  FLU_FV_22_${LOCI} $LOCI &
    bash make_reference_fake.sh /Nancy/data/input/ImmunoSeq/AbVitro/flu_time_course/FV/single_reads_with_UMIs/FV_23.fastq  FLU_FV_23_${LOCI} $LOCI &
    bash make_reference_fake.sh /Nancy/data/input/ImmunoSeq/AbVitro/flu_time_course/FV/single_reads_with_UMIs/FV_27.fastq  FLU_FV_27_${LOCI} $LOCI &
    bash make_reference_fake.sh /Engels/data/input/ImmunoSeq/Yale_UMI_Error_Correction/AAYHL_MG91M/AAYHL_MG91M-R12_primers-pass_pair-pass.fastq MG91M_${LOCI} $LOCI &
    bash make_reference_fake.sh /Engels/data/input/ImmunoSeq/Yale_UMI_Error_Correction/AAYHL_HD09M/AAYHL_HD09M-R12_primers-pass_pair-pass.fastq HD09M_${LOCI} $LOCI &
    wait
done

for i in 1 2 3 4 5 6 7 8 9
do
    bash make_reference_fake.sh /Nancy/data/input/ImmunoSeq/ibh_datasets/age_datasets/raw_reads/merged_reads/age_ig_s${i}_R12.fastq AGE${i} &
done
wait
