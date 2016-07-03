#!/bin/bash


INPUT=$1
OUTPUT=$2

IGREC_DIR=/home/ashlemov/Git/ig_repertoire_constructor

BCLEANER=${IGREC_DIR}/src/extra/ig_quast_tool/barcode_cleaner_nojoin.py
CFINDER="${IGREC_DIR}/build/release/bin/ig_consensus_finder -H"
CIC="${IGREC_DIR}/py/ig_compress_equal_clusters.py "

LOCI=all

${IGREC_DIR}/igrec.py -s ${INPUT} -o ${OUTPUT}/igrec_for_align/ --loci=${LOCI} --create-triv-dec -t4 --tau=1

mv ${OUTPUT}/igrec_for_align/vj_finder/cleaned_reads.fa ${OUTPUT}/cleaned_reads.fa
gzip ${OUTPUT}/cleaned_reads.fa -f
rm -fr ${OUTPUT}/igrec_for_align

${BCLEANER} ${OUTPUT}/cleaned_reads.fa.gz ${OUTPUT}/input1.fa.gz -r ${OUTPUT}/repertoire1.rcm --tau=0 -d100500 --distance-plot=${OUTPUT}/dist
${BCLEANER} ${OUTPUT}/cleaned_reads.fa.gz ${OUTPUT}/input2.fa.gz -r ${OUTPUT}/repertoire2.rcm --tau=2 -d100500
${BCLEANER} ${OUTPUT}/cleaned_reads.fa.gz ${OUTPUT}/input3.fa.gz -r ${OUTPUT}/repertoire3.rcm --tau=2 -d10

for i in 1 2 3
do
    ${CFINDER} -R ${OUTPUT}/repertoire${i}.rcm -i ${OUTPUT}/input${i}.fa.gz -o ${OUTPUT}/repertoire${i}.fa.gz
    ${CIC} -r ${OUTPUT}/repertoire${i}.rcm ${OUTPUT}/repertoire${i}.fa.gz ${OUTPUT}/repertoire${i}.fa.gz
done
