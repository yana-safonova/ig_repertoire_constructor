#!/bin/bash

# MG91M=/Engels/data/input/ImmunoSeq/Yale_UMI_Error_Correction/AAYHL_MG91M/AAYHL_MG91M-R12_primers-pass_pair-pass.fastq
HD09M=/Engels/data/input/ImmunoSeq/Yale_UMI_Error_Correction/AAYHL_HD09M/AAYHL_HD09M-R12_primers-pass_pair-pass.fastq

MUSCLE="/home/sbankevich/.local/bin/muscle"
USEARCH="/home/sbankevich/.local/bin/usearch"

set -e
ClusterSets.py -s $1 -f BARCODE -k CLUSTER --exec ${USEARCH} --nproc=16 --outname R12

ParseHeaders.py copy -s R12_cluster-pass.fasta -f BARCODE -k CLUSTER --act cat


BARCODE_NAME="CLUSTER"

BuildConsensus.py -s R12_cluster-pass_reheader.fasta --bf ${BARCODE_NAME} --pf PRIMER \
    --prcons 0.6 --maxerror 0.1 --maxgap 0.5 --outname MS12_R12 --log BC12.log --nproc=16

ParseHeaders.py collapse -s MS12_R12_consensus-pass.fasta -f CONSCOUNT --act min

CollapseSeq.py -s MS12_R12_consensus-pass_reheader.fasta -n 20 --inner --uf PRCONS \
    --cf CONSCOUNT --act sum --outname MS12



SplitSeq.py group -s MS12_collapse-unique.fasta -f CONSCOUNT --num 2 --outname MS12
SplitSeq.py group -s MS12_collapse-unique.fasta -f CONSCOUNT --num 5 --outname MS12
SplitSeq.py group -s MS12_collapse-unique.fasta -f DUPCOUNT --num 2 --outname MS12_BARCODES

ParseHeaders.py table -s MS12_atleast-5.fasta -f ID PRCONS CONSCOUNT DUPCOUNT

