#!/bin/bash

ref_AGE3="/home/ashlemov/Git/ig_repertoire_constructor/src/extra/ig_quast_tool/AGE3/repertoire3.fa.gz"
ref_AGE7="/home/ashlemov/Git/ig_repertoire_constructor/src/extra/ig_quast_tool/AGE7/repertoire3.fa.gz"
ref_SIM="/home/ashlemov/Git/ig_repertoire_constructor/various_error_rate/errate_0.0000/ideal_final_repertoire.fa.gz"
ref_FLU="/home/ashlemov/Git/ig_repertoire_constructor/var_err_rate_real/errate_0.0000/ideal_final_repertoire.fa.gz"

for file in "${ref_AGE3}" "${ref_AGE7}" "${ref_SIM}" "${ref_FLU}"
do
    ./repertoire_stats.py $file
done
    

./repertoire_stats.py "/home/ashlemov/Git/ig_repertoire_constructor/src/extra/ig_quast_tool/FLU_FV_21_IGL/repertoire3.fa.gz"
./repertoire_stats.py "/home/ashlemov/Git/ig_repertoire_constructor/src/extra/ig_quast_tool/FLU_FV_21_IG/repertoire3.fa.gz"
./repertoire_stats.py "/home/ashlemov/Git/ig_repertoire_constructor/src/extra/ig_quast_tool/FLU_FV_21_IG/repertoire1.fa.gz"
