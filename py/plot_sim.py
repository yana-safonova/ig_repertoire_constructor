#!/usr/bin/env python2

from plot_test import *

if __name__ == "__main__":
    rocs(igrec_dir + "/various_error_rate/errate_2.0000",
         tools=["igrec", "mixcr", "supernode"],
         labels=["IgReC", "MiXCR", "pRESTO"],
         title="SIMULATED, error rate = 2, with answer",
         out="sensitivity_precision_plot_er2")
    rocs(igrec_dir + "/various_error_rate/errate_1.0000",
         tools=["igrec", "mixcr", "supernode"],
         labels=["IgReC", "MiXCR", "pRESTO"],
         title="SIMULATED, error rate = 1, with answer",
         out="sensitivity_precision_plot_er1")
    rocs(igrec_dir + "/various_error_rate/errate_0.5000",
         tools=["igrec", "mixcr", "supernode"],
         labels=["IgReC", "MiXCR", "pRESTO"],
         title="SIMULATED, error rate = 0.5, with answer",
         out="sensitivity_precision_plot_er05")

    rocs(igrec_dir + "/various_error_rate/errate_2.0000_woans",
         tools=["igrec"],
         labels=["IgReC"],
         title="SIMULATED, error rate = 2, without answer",
         out="sensitivity_precision_plot_er2_woanswer")
    rocs(igrec_dir + "/various_error_rate/errate_1.0000_woans",
         tools=["igrec"],
         labels=["IgReC"],
         title="SIMULATED, error rate = 1, without answer",
         out="sensitivity_precision_plot_er1_woanswer")
    rocs(igrec_dir + "/various_error_rate/errate_0.5000_woans",
         tools=["igrec"],
         labels=["IgReC"],
         title="SIMULATED, error rate = 0.5, without answer",
         out="sensitivity_precision_plot_er05_woanswer")

    rocs(igrec_dir + "/v_e_r_r/errate_2.0000",
         tools=["igrec", "mixcr", "supernode"],
         labels=["IgReC", "MiXCR", "pRESTO"],
         title="SYNTHETIC, error rate = 2, with answer",
         out="sensitivity_precision_plot_er2_real")
    rocs(igrec_dir + "/v_e_r_r/errate_1.0000",
         tools=["igrec", "mixcr", "supernode"],
         labels=["IgReC", "MiXCR", "pRESTO"],
         title="SYNTHETIC, error rate = 1, with answer",
         out="sensitivity_precision_plot_er1_real")
    rocs(igrec_dir + "/v_e_r_r/errate_0.5000",
         tools=["igrec", "mixcr", "supernode"],
         labels=["IgReC", "MiXCR", "pRESTO"],
         title="SYNTHETIC, error rate = 0.5, with answer",
         out="sensitivity_precision_plot_er05_real")

    # for i, toolname in enumerate(["igrec", "mixcr", "presto"]):
    #     plot_various_error_rate(igrec_dir + "/various_error_rate", what="sensitivity", out="sensitivity_%s" % toolname,
    #                             which=[i],
    #                             kinds=["igrec_tau3_f03", "mixcr", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"],
    #                             legend=False)
    #     plot_various_error_rate(igrec_dir + "/various_error_rate", what="precision", out="precision_%s" % toolname,
    #                             which=[i],
    #                             kinds=["igrec_tau3_f03", "mixcr", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"],
    #                             legend=False)

    plot_various_error_rate(igrec_dir + "/various_error_rate", what="sensitivity", out="sensitivity",
                            title="SIMULATED, sensitivity, with answer",
                            kinds=["igrec", "mixcr", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])
    plot_various_error_rate(igrec_dir + "/various_error_rate", what="precision", out="precision",
                            title="SIMULATED, precision, with answer",
                            kinds=["igrec", "mixcr", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])
    plot_various_error_rate(igrec_dir + "/various_error_rate", out="sensitivity_plus_precision", what="sum",
                            title="SIMULATED, sensitivity + precision, with answer",
                            kinds=["igrec", "mixcr", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])
    plot_various_error_rate(igrec_dir + "/various_error_rate", out="optimal_minsize", what="minsize",
                            title="SIMULATED, optimal min size, with answer",
                            kinds=["igrec", "mixcr", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])


    plot_two_sums(igrec_dir + "/various_error_rate", out="two_sum",
                  title="SIMULATED, sens + prec, with/without answer",
                  kind="igrec", label="IgReC")


    # plot_various_error_rate(igrec_dir + "/various_error_rate", what="sensitivity", out="sensitivity_woans", woans=True,
    #                         title="SIMULATED, sensitivity, without answer",
    #                         kinds=["igrec_tau3", "mixcr", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])
    # plot_various_error_rate(igrec_dir + "/various_error_rate", what="precision", out="precision_woans", woans=True,
    #                         title="SIMULATED, precision, without answer",
    #                         kinds=["igrec_tau3", "mixcr", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])

    plot_various_error_rate(igrec_dir + "/v_e_r_r", what="sensitivity", out="sensitivity_flu",
                            title="SYNTHETIC FLU, sensitivity",
                            kinds=["igrec", "mixcr", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])
    plot_various_error_rate(igrec_dir + "/v_e_r_r", what="precision", out="precision_flu",
                            title="SYNTHETIC FLU, precision",
                            kinds=["igrec", "mixcr", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])
    plot_various_error_rate(igrec_dir + "/v_e_r_r", out="sensitivity_plus_precision_flu", what="sum",
                            title="SYNTHETIC FLU, sensitivity + precision",
                            kinds=["igrec", "mixcr", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])

    # plot_various_error_rate(igrec_dir + "/various_error_rate",
    #                         kind1="igrec_tau3_f03", kind2="mixcr",
    #                         what="precision", out="precision_tau3")
    # plot_various_error_rate(igrec_dir + "/various_error_rate",
    #                         kind1="igrec_tau3_f03", kind2="mixcr",
    #                         what="sensitivity", out="sensitivity_tau3")
