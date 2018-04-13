#!/bin/bash

set -euxo pipefail

# FROM https://stackoverflow.com/questions/1537956/bash-limit-the-number-of-concurrent-jobs
function max_bg_procs {
    if [[ $# -eq 0 ]] ; then
            echo "Usage: max_bg_procs NUM_PROCS.  Will wait until the number of background (&)"
            echo "           bash processes (as determined by 'jobs -pr') falls below NUM_PROCS"
            return
    fi
    local max_number=$((0 + ${1:-0}))
    set +x
    while true; do
            local current_number=$(jobs -pr | wc -l)
            if [[ $current_number -lt $max_number ]]; then
                    break
            fi
            sleep 1
    done
    set -x
}


# for LOCI in IG IGH IGL IGK
# do
#     bash make_reference.sh /Nancy/data/input/ImmunoSeq/AbVitro/flu_time_course/FV/single_reads_with_UMIs/FV_21.fastq  FLU_FV_21_${LOCI} $LOCI &
#     bash make_reference.sh /Nancy/data/input/ImmunoSeq/AbVitro/flu_time_course/FV/single_reads_with_UMIs/FV_22.fastq  FLU_FV_22_${LOCI} $LOCI &
#     bash make_reference.sh /Nancy/data/input/ImmunoSeq/AbVitro/flu_time_course/FV/single_reads_with_UMIs/FV_23.fastq  FLU_FV_23_${LOCI} $LOCI &
#     bash make_reference.sh /Nancy/data/input/ImmunoSeq/AbVitro/flu_time_course/FV/single_reads_with_UMIs/FV_27.fastq  FLU_FV_27_${LOCI} $LOCI &
#     bash make_reference.sh /Engels/data/input/ImmunoSeq/Yale_UMI_Error_Correction/AAYHL_MG91M/AAYHL_MG91M-R12_primers-pass_pair-pass.fastq MG91M_${LOCI} $LOCI &
#     bash make_reference.sh /Engels/data/input/ImmunoSeq/Yale_UMI_Error_Correction/AAYHL_HD09M/AAYHL_HD09M-R12_primers-pass_pair-pass.fastq HD09M_${LOCI} $LOCI &
#     wait
# done
#
for i in 1 2 3 4 5 6 7 8 9
do
    max_bg_procs 4
    bash make_reference.sh /Nancy/data/input/ImmunoSeq/ibh_datasets/age_datasets/raw_reads/merged_reads/age_ig_s${i}_R12.fastq AGE${i} &
done
wait
#
#

# for i in 1 2
# do
#     for j in 1 2
#     do
#         bash make_reference.sh ~/igrec_revisited/donor1_S1.fasta.gz donor_${i}_S${j} &
#     done
# done
# wait


# for id in `cat ids.txt`
# do
#   max_bg_procs 4
#   bash make_reference.sh ~/igrec_revisited/Stern/finished/barcoded_${id}.fastq.gz Stern_${id} IG &
# done
# wait
