#!/bin/bash


INPUT=$1
OUTPUT=$2
OUTPUT_READS="${OUTPUT}/input_reads.fa.gz"
OUTPUT_CONSENSUS="${OUTPUT}/repertoire.fa.gz"
OUTPUT_RCM="${OUTPUT}/repertoire.rcm"
OUTPUT_CONSENSUS_SHORT="${OUTPUT}/repertoire_short.fa.gz"
OUTPUT_RCM_SHORT="${OUTPUT}/repertoire_short.rcm"

IGREC_DIR=/home/ashlemov/Git/ig_repertoire_constructor


${IGREC_DIR}/igrec.py -s ${INPUT} -o ${OUTPUT}/igrec/ --loci=all --create-triv-dec

${IGREC_DIR}/src/extra/ig_quast_tool/barcode_cleaner_nojoin.py ${OUTPUT}/igrec/vj_finder/cleaned_reads.fa ${OUTPUT_READS} -r ${OUTPUT_RCM}
cp ${OUTPUT_RCM} ${OUTPUT_RCM_SHORT}

${IGREC_DIR}/build/release/bin/ig_consensus_finder -H -R ${OUTPUT_RCM} -i ${OUTPUT_READS} -o ${OUTPUT_CONSENSUS} -l 0
${IGREC_DIR}/py/ig_compress_equal_clusters.py ${OUTPUT_CONSENSUS} ${OUTPUT_CONSENSUS} --rcm ${OUTPUT_RCM} --output-rcm ${OUTPUT_RCM}

${IGREC_DIR}/build/release/bin/ig_consensus_finder -H -R ${OUTPUT_RCM_SHORT} -i ${OUTPUT_READS} -o ${OUTPUT_CONSENSUS_SHORT}
${IGREC_DIR}/py/ig_compress_equal_clusters.py ${OUTPUT_CONSENSUS_SHORT} ${OUTPUT_CONSENSUS_SHORT} --rcm ${OUTPUT_RCM_SHORT} --output-rcm ${OUTPUT_RCM_SHORT}


${IGREC_DIR}/igrec.py -s ${OUTPUT_READS} -o ${OUTPUT}/igrec_cleaned/ --loci=all
