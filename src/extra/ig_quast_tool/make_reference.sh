#!/bin/bash


INPUT=$1
OUTPUT=$2
OUTPUT_READS="${OUTPUT}/input_reads.fa.gz"
OUTPUT_CONSENSUS="${OUTPUT}/repertoire.fa.gz"
OUTPUT_RCM="${OUTPUT}/repertoire.rcm"

IGREC_DIR=/home/ashlemov/Git/ig_repertoire_constructor


# ${IGREC_DIR}/igrec.py -s ${INPUT} -o ${OUTPUT}/igrec/ --loci=all --create-triv-dec
#
# ${IGREC_DIR}/src/extra/ig_quast_tool/barcode_cleaner_nojoin.py ${OUTPUT}/igrec/vj_finder/cleaned_reads.fa ${OUTPUT_READS} -r ${OUTPUT_RCM}
${IGREC_DIR}/src/extra/ig_quast_tool/barcode_cleaner_nojoin.py ${INPUT} ${OUTPUT_READS} -r ${OUTPUT_RCM}

${IGREC_DIR}/build/release/bin/ig_consensus_finder -H -R ${OUTPUT_RCM} -i ${OUTPUT_READS} -o ${OUTPUT_CONSENSUS}
${IGREC_DIR}/py/ig_compress_equal_clusters.py ${OUTPUT_CONSENSUS} ${OUTPUT_CONSENSUS} --rcm ${OUTPUT_RCM} --output-rcm ${OUTPUT_RCM}


rm -fr ${OUTPUT}/igrec
${IGREC_DIR}/igrec.py -s ${OUTPUT_READS} -o ${OUTPUT}/igrec/ --loci=all



rm -fr ${OUTPUT}/igrec/vj_finder

rm ${OUTPUT}/igrec/super_reads.fa
rm ${OUTPUT}/igrec/final_repertoire_large.fa

gzip ${OUTPUT}/igrec/final_repertoire.fa


${IGREC_DIR}/igrec.py -s ${OUTPUT_READS} -o ${OUTPUT}/igrec_trivial/ --loci=all --create-triv-dec

rm -fr ${OUTPUT}/igrec_trivial/vj_finder
rm ${OUTPUT}/igrec_trivial/super_reads.fa
rm ${OUTPUT}/igrec_trivial/final_repertoire_large.fa
gzip ${OUTPUT}/igrec_trivial/final_repertoire.fa
