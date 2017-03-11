#!/usr/bin/env python2

from plot_test import *


def plot_sens_prec_umi(base_results_dir, pcr_error_rates = (0.0006, 0.0025, 0.006), error_rates = (0.5, 2, 5)):
    lambdas = pcr_error_rates
    lambdas_str = ["low", "medium", "high"]
    lambdas_str.extend('?' * (len(pcr_error_rates) - 3))
    for lam, lam_str, error_rate in zip(lambdas, lambdas_str, error_rates):
        rocs("%s/pcr_%g_super_100000_umi_15/" % (base_results_dir, lam),
            tools=["quast_barigrec", "quast_igrec", "quast_presto", "quast_migec"],
            labels=["BarcodedIgReC", "IgReC", "pRESTO", "MiGEC + MiXCR"],
            # title="test_plot_fot serg",
            title="PCR error rate = %0.4f, read error rate = %0.1f" % (lam, error_rate),
             out="%s/plots/barigrec_%f" % (base_results_dir, lam),
             add_aimquast_to_path=False)

def plot_sens_prec_umi_tau(base_results_dir, min_tau = 0, max_tau = 3, pcr_error_rates = (0.006, 0.0006, 0.0025, 0.0012, 0.0018, 0.0030, 0.0036)):
    taus = range(min_tau, max_tau + 1)
    for lam in pcr_error_rates:
        rocs("%s/pcr_%1.6f_super_100000_umi_15/" % (base_results_dir, lam),
            tools=["tau_%d" % tau for tau in taus],
            labels=["tau_%d" % tau for tau in taus],
            # title="test_plot_fot serg",
            # title="PCR error rate = %0.4f, read error rate = %0.1f" % (lam, error_rate),
             out="%s/plots/barigrec_%f_tau" % (base_results_dir, lam),
             add_aimquast_to_path=False)


if __name__ == "__main__":
    # plot_sens_prec_umi(sys.argv[1])
    plot_various_error_rate_serg("/Marx/serg/data/ig_simulator/new_error_rates/",
                                 what="sensitivity", out="sensitivity_fig",
                                 title="Sensitivity plot\n(SIMULATED BARCODED dataset)",
                                 kinds=["barigrec", "igrec", "presto", "migec"], labels=["BarcodedIgReC", "IgReC", "pRESTO", "MiGEC + MiXCR"])
    plot_various_error_rate_serg("/Marx/serg/data/ig_simulator/new_error_rates/",
                                 what="precision", out="precision_fig",
                                 title="Precision plot\n(SIMULATED BARCODED dataset)",
                                 kinds=["barigrec", "igrec", "presto", "migec"], labels=["BarcodedIgReC", "IgReC", "pRESTO", "MiGEC + MiXCR"])
    plot_various_error_rate_serg("/Marx/serg/data/ig_simulator/new_error_rates/",
                                 what="sum", out="sens_prec_fig",
                                 title="Sensitivity * precision plot\n(SIMULATED BARCODED dataset)",
                                 legend_loc=2,
                                 kinds=["barigrec", "igrec", "presto", "migec"], labels=["BarcodedIgReC", "IgReC", "pRESTO", "MiGEC + MiXCR"])
    plot_various_error_rate_serg("/Marx/serg/data/ig_simulator/new_error_rates/",
                                 what="minsize", out="min_size_fig",
                                 title="Size threshold plot\n(SIMULATED BARCODED dataset)",
                                 # title="Size threshold plot (SIMULATED BARCODED dataset, 2.0 errors per read)",
                                 legend_loc=2,
                                 kinds=["barigrec", "igrec", "presto", "migec"], labels=["BarcodedIgReC", "IgReC", "pRESTO", "MiGEC + MiXCR"])
