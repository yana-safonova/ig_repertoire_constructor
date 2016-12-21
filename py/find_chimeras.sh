#!/bin/bash

set -e

INPUT_FILE=$1
OUTDIR=$2

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

mkdir -p $OUTDIR

${DIR}/convert_repertoire_for_uchime.py ${INPUT_FILE} ${OUTDIR}/converted.fa

usearch -uchime2_denovo ${OUTDIR}/converted.fa -uchimeout ${OUTDIR}/out.txt -chimeras ${OUTDIR}/ch.fa -nonchimeras ${OUTDIR}/nonch.fa
