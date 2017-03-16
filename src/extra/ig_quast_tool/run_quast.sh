#!/bin/bash


INPUT=$1
INPUT_READS="${INPUT}/input_reads.fa.gz"
REFERENCE_CONSENSUS="${INPUT}/repertoire.fa.gz"
REFERENCE_RCM="${INPUT}/repertoire.rcm"
REFERENCE_CONSENSUS_SHORT="${INPUT}/repertoire_short.fa.gz"
REFERENCE_RCM_SHORT="${INPUT}/repertoire_short.rcm"
IGREC_CONSENSUS="${INPUT}/igrec_cleaned/final_repertoire.fa.gz"
IGREC_RCM="${INPUT}/igrec_cleaned/final_repertoire.rcm"
IGREC_CONSENSUS_TRIVIAL="${INPUT}/igrec_trivial/final_repertoire.fa.gz"
IGREC_RCM_TRIVIAL="${INPUT}/igrec_trivial/final_repertoire.rcm"


IGREC_DIR=/home/ashlemov/Git/ig_repertoire_constructor

OUTPUT="${INPUT}/quast"

${IGREC_DIR}/igquast.py -s ${INPUT_READS} -c ${IGREC_CONSENSUS} -C ${IGREC_RCM} -r ${REFERENCE_CONSENSUS} -R ${REFERENCE_RCM} \
    -o ${OUTPUT} --json ${OUTPUT}/report.json --text ${OUTPUT}/report.txt \
    -F='svg,png,pdf'


OUTPUT_SHORT="${INPUT}/quast_short"

# ${IGREC_DIR}/igquast.py -s ${INPUT_READS} -c ${IGREC_CONSENSUS} -C ${IGREC_RCM} -r ${REFERENCE_CONSENSUS_SHORT} -R ${REFERENCE_RCM_SHORT} \
#     -o ${OUTPUT_SHORT} --json ${OUTPUT_SHORT}/report.json --text ${OUTPUT_SHORT}/report.txt \
#     -F='svg,png,pdf'


OUTPUT_SHORT_TRIVIAL="${INPUT}/quast_short_trivial"

${IGREC_DIR}/igquast.py -s ${INPUT_READS} -c ${IGREC_CONSENSUS_TRIVIAL} -C ${IGREC_RCM_TRIVIAL} -r ${REFERENCE_CONSENSUS_SHORT} -R ${REFERENCE_RCM_SHORT} \
    -o ${OUTPUT_SHORT_TRIVIAL} --json ${OUTPUT_SHORT_TRIVIAL}/report.json --text ${OUTPUT_SHORT_TRIVIAL}/report.txt \
    -F='svg,png,pdf'


