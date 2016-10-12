#!/usr/bin/env python2

from plot_test import *

if __name__ == "__main__":
    lambdas = [0.0002, 0.0006, 0.002]
    lambdas_str = ["low", "medium", "high"]
    for lam, lam_str in zip(lambdas, lambdas_str):
        rocs("/Marx/serg/data/ig_simulator/pcr_%1.6f_super_100000_umi_15/" % lam,
            tools=["quast", "quast_igrec", "quast_presto", "quast_migec"],
            labels=["barcoded IgReC", "IgReC", "pRESTO", "MiGEC"],
            # title="test_plot_fot serg",
            title="%s PCR error rate" % lam_str,
            out="/Marx/SergAndAlex/plots/barigrec_%f" % lam,
            add_aimquast_to_path=False)
