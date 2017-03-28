#!/bin/bash


INPUT=$1
OUTPUT=$2
set -e

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

IGREC_DIR=${DIR}/../../../.

BCLEANER=${IGREC_DIR}/src/extra/ig_quast_tool/barcode_cleaner_nojoin.py
CFINDER="${IGREC_DIR}/build/release/bin/ig_component_splitter"
CIC="${IGREC_DIR}/py/ig_compress_equal_clusters.py "

LOCI=${3:-all}

${IGREC_DIR}/igrec.py -s ${INPUT} -o ${OUTPUT}/igrec_for_align/ --loci=${LOCI} --create-triv-dec -t4 --tau=1

mv ${OUTPUT}/igrec_for_align/vj_finder/cleaned_reads.fa ${OUTPUT}/cleaned_reads.fa
gzip ${OUTPUT}/cleaned_reads.fa -f
rm -fr ${OUTPUT}/igrec_for_align

${BCLEANER} ${OUTPUT}/cleaned_reads.fa.gz ${OUTPUT}/input1.fa.gz -r ${OUTPUT}/repertoire1.rcm --tau=0 -d100500 --distance-plot=${OUTPUT}/dist --lengths ${OUTPUT}/lengths.txt
${BCLEANER} ${OUTPUT}/cleaned_reads.fa.gz ${OUTPUT}/input2.fa.gz -r ${OUTPUT}/repertoire2.rcm --tau=2 -d100500
${BCLEANER} ${OUTPUT}/cleaned_reads.fa.gz ${OUTPUT}/input3.fa.gz -r ${OUTPUT}/repertoire3.rcm --tau=2 -d10



for i in 1 2 3
do
    ${CFINDER} -R ${OUTPUT}/repertoire${i}.rcm -i ${OUTPUT}/input${i}.fa.gz -o ${OUTPUT}/repertoire${i}.fa.gz
    ${CIC} -r ${OUTPUT}/repertoire${i}.rcm ${OUTPUT}/repertoire${i}.fa.gz ${OUTPUT}/repertoire${i}.fa.gz -S ${OUTPUT}/barcode_stats${i}.tsv
done


# TODO reconstruct input4
# cp ${OUTPUT}/input3.fa.gz ${OUTPUT}/input4.fa.gz
# cp ${OUTPUT}/repertoire3.rcm ${OUTPUT}/repertoire4.rcm
# cp ${OUTPUT}/repertoire3.fa.gz ${OUTPUT}/repertoire4.fa.gz
# ${CIC} -r ${OUTPUT}/repertoire4.rcm ${OUTPUT}/repertoire4.fa.gz ${OUTPUT}/repertoire4.fa.gz -b 2


JITTER=${IGREC_DIR}/py/jit_file.py
for lambda in "0.25" "0.5" "1" "2"
do
    for i in 1 2 3 4
    do
        ${JITTER} ${OUTPUT}/input${i}.fa.gz ${OUTPUT}/input${i}_jit${lambda}.fa.gz --error-rate=${lambda} &
    done
    wait
done
