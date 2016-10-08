import os
import shutil
import sys

from string import join

from os import getcwd


def run_sim_pipeline(data_path, pcr_error_rate, supernode_threshold, barcode_length, threads):
    home = "~/work" if os.path.exists("~/work") else "~"
    igrec = home + "/igrec"
    args_seq = [
        ["%s/build/release/bin/simulate_barcoded" % igrec,
         "--input-file %s/final_repertoire.fasta" % data_path,
         "--output-file %s/amplified.fasta" % data_path,
         "--umi-length %d" % barcode_length,
         "--pcr-error1 %f" % pcr_error_rate,
         "--pcr-error2 %f" % pcr_error_rate,
         "--compressed-path %s/final_repertoire_comp.fasta" % data_path
         ],
        ["cd %s &&" % igrec,
         "%s/build/release/bin/vj_finder" % igrec,
         "--input-file %s/final_repertoire_comp.fasta" % data_path,
         "--output-dir %s/vjf" % data_path,
         "--loci IG",
         "--threads %d" % threads
         ],

        ["python %s/igrec_umi.py" % igrec,
         "-s %s/amplified.fasta" % data_path,
         "--output %s/igrec_umi" % data_path,
         "--loci IG",
         "--threads %d" % threads,
         "--igrec_tau 2",
         "--min-super-read-size %d" % supernode_threshold,
         "--no-compilation",
         "--detect-chimeras",
         "--clustering-thr 20"
         ],
        ["python %s/py/drop_ns.py" % igrec,
         "%s/igrec_umi/final_repertoire/final_repertoire.fa",
         "%s/igrec_umi/final_repertoire.fa",
         "%s/igrec_umi/final_repertoire/final_repertoire.rcm",
         "%s/igrec_umi/final_repertoire.rcm",
         ],
        ["python %s/aimquast.py" % igrec,
         "-s %s/amplified.fasta" % data_path,
         "-c %s/igrec_umi/final_repertoire.fa" % data_path,
         "-C %s/igrec_umi/final_repertoire.rcm" % data_path,
         "-r %s/vjf/cleaned_reads.fa" % data_path,
         "-o %s/quast" % data_path,
         "--reference-free",
         "--rcm-based"
         ],
        # ["python %s/aimquast.py" % igrec,
        #  "-s %s/amplified.fasta" % data_path,
        #  "-c %s/igrec_umi/final_repertoire/final_repertoire.fa" % data_path,
        #  "-C %s/igrec_umi/final_repertoire/final_repertoire.rcm" % data_path,
        #  "-r %s/vjf/cleaned_reads.fa" % data_path,
        #  "-o %s/quast" % data_path,
        #  "--reference-free",
        #  "--rcm-based"
        #  ],

        ["python %s/igrec.py" % igrec,
         "-s %s/amplified.fasta" % data_path,
         "-o %s/igrec" % data_path,
         "--threads %d" % threads,
         "--loci IGH",
         "--debug"
         ],
        ["python %s/py/drop_ns.py" % igrec,
         "%s/igrec/final_repertoire.fa",
         "%s/igrec/final_repertoire_non.fa",
         "%s/igrec/final_repertoire.rcm",
         "%s/igrec/final_repertoire_non.rcm",
         ],
        ["python %s/aimquast.py" %igrec,
         "-s %s/amplified.fasta" % data_path,
         "-c %s/igrec/final_repertoire_non.fa" % data_path,
         "-C %s/igrec/final_repertoire_non.rcm" % data_path,
         "-r %s/vjf/cleaned_reads.fa" % data_path,
         "-o %s/quast_igrec" % data_path,
         "--reference-free",
         "--rcm-based"
         ],
        # ["python %s/aimquast.py" %igrec,
        #  "-s %s/amplified.fasta" % data_path,
        #  "-c %s/igrec/final_repertoire.fa" % data_path,
        #  "-C %s/igrec/final_repertoire.rcm" % data_path,
        #  "-r %s/vjf/cleaned_reads.fa" % data_path,
        #  "-o %s/quast_igrec" % data_path,
        #  "--reference-free",
        #  "--rcm-based"
        #  ],

        ["mkdir -p %s/presto &&" % data_path,
         "python /Marx/ashlemov/Git/ig_repertoire_constructor/py/convertAGE2PRESTO.py",
         "%s/amplified.fasta" % data_path,
         "%s/presto/amplified_for_presto.fasta" % data_path
        ],
        ["cd %s/presto &&" % data_path,
         "../../run_simple.sh",
         "amplified_for_presto.fasta"
         ],
        ["python %s/py/convert_presto_to_quast.py" % igrec,
         "-r %s/presto/MS12_collapse-unique.fasta" % data_path,
         "-o %s/presto/presto.fasta" % data_path
         ],
        ["python %s/py/drop_ns.py" % igrec,
         "%s/presto/presto.fasta",
         "%s/presto/presto_non.fasta"
         ],
        ["python %s/aimquast.py" %igrec,
         "-s %s/amplified.fasta" % data_path,
         "-c %s/presto/presto_non.fasta" % data_path,
         "-r %s/vjf/cleaned_reads.fa" % data_path,
         "-o %s/quast_presto" % data_path,
         "--reference-free",
         "--rcm-based"
         ],
        # ["python %s/aimquast.py" %igrec,
        #  "-s %s/amplified.fasta" % data_path,
        #  "-c %s/presto/presto.fasta" % data_path,
        #  "-r %s/vjf/cleaned_reads.fa" % data_path,
        #  "-o %s/quast_presto" % data_path,
        #  "--reference-free",
        #  "--rcm-based"
        #  ]
    ]

    commands = [join(args, ' ') for args in args_seq]
    for command in commands:
        print "Running %s" % command
        exit_status = os.system(command)
        print "Returned", exit_status
        if exit_status != 0:
            exit(exit_status)


def main():
    for pcr_error_rate in [0.002, 0.0006, 0.0002]:
        for supernode_threshold in [5, 10, 100000]:
            for barcode_length in [9, 12, 15]:
                data_path = "%s/pcr_%f_super_%d_umi_%d" % (os.getcwd(), pcr_error_rate, supernode_threshold, barcode_length)
                if not os.path.exists(data_path):
                    os.makedirs(data_path)
                shutil.copyfile("final_repertoire.fasta", "%s/final_repertoire.fasta" % data_path)
                run_sim_pipeline(data_path, pcr_error_rate, supernode_threshold, barcode_length, int(sys.argv[1]))


if __name__ == '__main__':
    main()
