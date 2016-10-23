#!/bin/bash

set -e

./find_chimeras.sh ../src/extra/ig_quast_tool/AGE3/repertoire3.fa.gz AGE3_chimeras &
./find_chimeras.sh ../src/extra/ig_quast_tool/AGE7/repertoire3.fa.gz AGE7_chimeras &
./find_chimeras.sh ../src/extra/ig_quast_tool/FLU_FV_21_IG/repertoire3.fa.gz FLU_FV_21_chimeras &
./find_chimeras.sh ../src/extra/ig_quast_tool/MG91M_IG/repertoire3.fa.gz MG91M_chimeras &
./find_chimeras.sh ../src/extra/ig_quast_tool/HD09M_IG//repertoire3.fa.gz HD09M_chimeras &
wait

