#!/usr/bin/env python2

from plot_test import *


def plot_sens_prec_umi(base_results_dir):
    lambdas = [0.0006, 0.0025, 0.006]
    lambdas_str = ["low", "medium", "high"]
    for lam, lam_str in zip(lambdas, lambdas_str):
        rocs("%s/pcr_%1.6f_super_100000_umi_15/" % (base_results_dir, lam),
             tools=["quast", "quast_igrec", "quast_presto", "quast_migec"],
             labels=["barcoded IgReC", "IgReC", "pRESTO", "MiGEC"],
             # title="test_plot_fot serg",
             title="%s PCR error rate" % lam_str,
             out="%s/plots/barigrec_%f" % (sys.argv[1], lam),
             add_aimquast_to_path=False)


if __name__ == "__main__":
    plot_sens_prec_umi(sys.argv[1])
