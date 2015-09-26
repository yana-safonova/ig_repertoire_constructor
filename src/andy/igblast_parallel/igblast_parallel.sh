#!/bin/bash

IGBLASTDIR='/home/ashlemov/Yandex.Disk/Documents/lab/vers2/src/tools/igblast'
#directory of IMGT files
GERMLINE_DB='/home/ashlemov/Yandex.Disk/Documents/lab/vers2/src/tools/igblast/database'

INPUT_FILE=$1
OUTPUT_DIR=$2

if [ "${INPUT_FILE}" == "" ] || [ "${OUTPUT_DIR}" == "" ] || [ "${INPUT_FILE}" == "--help" ]; then
  echo "./igblast_parallel.sh <input_file> <output_file> <number_of_threads> <sequence format>"
  echo '      <number_of_threads> and <sequence format> are optional'
  echo '      for <number_of_threads> 8 is the default'
  echo '      For help write "./igblast_parallel.sh --help"' 
  exit 1
fi

OUTPUT_FILENAME=$(basename "$OUTPUT_DIR")
OUTPUT_DIR=$(dirname "$OUTPUT_DIR")

TEMP_DIR=`mktemp -d`

N_THREADS=$3
FORMAT=$4

export IGDATA=${IGBLASTDIR}

PATH_TO_SPLITTER='/home/ashlemov/Yandex.Disk/Documents/lab/vers2/src/andy/igblast_parallel/cut_merged_reads.py'

split="${PATH_TO_SPLITTER} -i $INPUT_FILE -o $TEMP_DIR"
#echo "N_Threads: ${N_THREADS}"

if [ "${N_THREADS}" != "" ]; then
  split="$split -n $N_THREADS"
fi

if [ "${FORMAT}" != "" ]; then
  split="$split -f $FORMAT"
fi

#echo "SPLIT: ${split}"
${split}

#ls $OUTPUT_DIR

for file in ${TEMP_DIR}/merged_reads_cut_???.fasta; do
  #echo ${file}
  ${IGBLASTDIR}/bin/igblastn -germline_db_V "${GERMLINE_DB}/human_gl_V" -germline_db_J "${GERMLINE_DB}/human_gl_J" -germline_db_D "${GERMLINE_DB}/human_gl_D" -domain_system imgt -query "${file}" -auxiliary_data "auxilary_file" -outfmt 7 -num_threads 1 -num_alignments_V 10 -out "${file}.blast" -show_translation &
done
wait

cat ${TEMP_DIR}/*.blast > "$OUTPUT_DIR/$OUTPUT_FILENAME"
rm -r $TEMP_DIR 

#igblastn -germline_db_V "${GERMLINE_DB}/human_gl_V" -germline_db_J "${GERMLINE_DB}/human_gl_J" -germline_db_D "${GERMLINE_DB}/human_gl_D" -domain_system imgt -query ${INPUT_FILE} -auxiliary_data "auxilary_file" -outfmt 5 -num_threads 2 -num_alignments_V 10 -out ${OUTPUT_FILE} -show_translation
