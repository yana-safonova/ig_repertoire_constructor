#!/bin/bash

get_script_dir () {
  #FROM: http://www.ostricher.com/2014/10/the-right-way-to-get-the-directory-of-a-bash-script/
  SOURCE="${BASH_SOURCE[0]}"
  # While $SOURCE is a symlink, resolve it
  while [ -h "$SOURCE" ]; do
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$( readlink "$SOURCE" )"
    # If $SOURCE was a relative symlink (so no "/" as prefix, need to resolve it relative to the symlink base directory
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
  done
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  echo "$DIR"
}

CUR_DIR="$(get_script_dir)"

IGREC_DIR=${CUR_DIR}/../../../

IGDATA_DIR="${CUR_DIR}/igblast"
IGBLASTDIR="${CUR_DIR}/igblast"
#directory of IMGT files

VJFINDER_GERMLINE_DB=${IGREC_DIR}/data/germline/human/IG
GERMLINE_DB=${IGDATA_DIR}/database

mkdir -p ${GERMLINE_DB}

(
  flock -n 9 || exit 1
  # ... commands executed under lock ...
  # Prepare database
  WORKDIR=`pwd`
  cd ${GERMLINE_DB}
  MAKEBLASTDB="makeblastdb"
  FORMATDB="formatdb"

  cp ${VJFINDER_GERMLINE_DB}/IGH*.fa ${GERMLINE_DB}

  for i in IGHD*.fa
  do
    $FORMATDB -o T -p F -i $i
  done

  for i in IGH{V,J}*.fa
  do
    $MAKEBLASTDB -parse_seqids -dbtype nucl -in $i
  done
  cd ${WORKDIR}
) 9>/var/lock/blastlock

INPUT_FILE=$1
OUTPUT_FILE=$2

export IGDATA=${IGDATA_DIR}

# ${IGBLASTDIR}/bin/igblastn \
# -auxiliary_data "${IGBLASTDIR}/auxiliary_file" \
# -germline_db_V "${GERMLINE_DB}/IGHV-allP.fa" \
# -germline_db_J "${GERMLINE_DB}/IGHJ-allP.fa" \
# -germline_db_D "${GERMLINE_DB}/IGHD-allP.fa" \
# -domain_system imgt -query ${INPUT_FILE} \
# -outfmt 7 -num_threads 1 -num_alignments_V 1  -num_alignments_D 0 -num_alignments_J 1 -organism human \
# -ig_seqtype Ig \
# -out ${OUTPUT_FILE}
${IGBLASTDIR}/bin/igblastn \
-auxiliary_data "${IGBLASTDIR}/auxiliary_file" \
-germline_db_V "${GERMLINE_DB}/IGHV.fa" \
-germline_db_J "${GERMLINE_DB}/IGHJ.fa" \
-germline_db_D "${GERMLINE_DB}/IGHD.fa" \
-domain_system imgt -query ${INPUT_FILE} \
-outfmt 7 -num_threads 1 -num_alignments_V 1  -num_alignments_D 0 -num_alignments_J 1 -organism human \
-ig_seqtype Ig \
-out ${OUTPUT_FILE}
