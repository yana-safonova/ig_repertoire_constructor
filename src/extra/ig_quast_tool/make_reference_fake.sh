#!/bin/bash


INPUT=$1
OUTPUT=$2

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

${BCLEANER} ${OUTPUT}/cleaned_reads.fa.gz ${OUTPUT}/input_fake.fa.gz -r ${OUTPUT}/repertoire_fake.rcm --tau=0 -d100500 --lengths ${OUTPUT}/lengths.txt.gz
rm ${OUTPUT}/input_fake.fa.gz
rm ${OUTPUT}/repertoire_fake.rcm
