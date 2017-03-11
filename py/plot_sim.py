#!/usr/bin/env python2

from plot_test import *

if __name__ == "__main__":

    # for i, toolname in enumerate(["igrec", "mixcr2", "presto"]):
    #     plot_various_error_rate(igrec_dir + "/SIMULATED", what="sensitivity", out="sensitivity_%s" % toolname,
    #                             which=[i],
    #                             kinds=["igrec_tau3_f03", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"],
    #                             legend=False)
    #     plot_various_error_rate(igrec_dir + "/SIMULATED", what="precision", out="precision_%s" % toolname,
    #                             which=[i],
    #                             kinds=["igrec_tau3_f03", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"],
    #                             legend=False)

    plot_various_error_rate(igrec_dir + "/SIMULATED", what="sensitivity", out="Fig_9a",
                            title="Sensitivity (SIMULATED SIMPLE)",
                            kinds=["igrec", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])
    plot_various_error_rate(igrec_dir + "/SIMULATED", what="precision", out="Fig_9b",
                            title="Precision (SIMULATED SIMPLE)",
                            kinds=["igrec", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])

    plot_various_error_rate(igrec_dir + "/SIMULATED", what="sensitivity", out="Fig_9a_prod",
                            prod_criterion=True,
                            title="Sensitivity (SIMULATED SIMPLE)",
                            kinds=["igrec", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])
    plot_various_error_rate(igrec_dir + "/SIMULATED", what="precision", out="Fig_9b_prod",
                            prod_criterion=True,
                            title="Precision (SIMULATED SIMPLE)",
                            kinds=["igrec", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])


    plot_various_error_rate(igrec_dir + "/SIMULATED", what="sensitivity", out="Fig_9a_with_IGRC",
                            title="Sensitivity (SIMULATED SIMPLE)",
                            kinds=["igrec", "mixcr2", "supernode", "ig_repertoire_constructor"], labels=["IgReC", "MiXCR", "pRESTO", "IgRepertoireConstructor"])
    plot_various_error_rate(igrec_dir + "/SIMULATED", what="precision", out="Fig_9b_with_IGRC",
                            title="Precision (SIMULATED SIMPLE)",
                            kinds=["igrec", "mixcr2", "supernode", "ig_repertoire_constructor"], labels=["IgReC", "MiXCR", "pRESTO", "IgRepertoireConstructor"])
    plot_various_error_rate(igrec_dir + "/SIMULATED", what="sum", out="Fig_9c_with_IGRC",
                            prod_criterion=True,
                            title="Sensitivity + precision (SIMULATED SIMPLE)",
                            kinds=["igrec", "mixcr2", "supernode", "ig_repertoire_constructor"], labels=["IgReC", "MiXCR", "pRESTO", "IgRepertoireConstructor"])

    plot_various_error_rate(igrec_dir + "/SIMULATED", what="sensitivity", out="Fig_9a_prod_with_IGRC",
                            prod_criterion=True,
                            title="Sensitivity (SIMULATED SIMPLE)",
                            kinds=["igrec", "mixcr2", "supernode", "ig_repertoire_constructor"], labels=["IgReC", "MiXCR", "pRESTO", "IgRepertoireConstructor"])
    plot_various_error_rate(igrec_dir + "/SIMULATED", what="precision", out="Fig_9b_prod_with_IGRC",
                            prod_criterion=True,
                            title="Precision (SIMULATED SIMPLE)",
                            kinds=["igrec", "mixcr2", "supernode", "ig_repertoire_constructor"], labels=["IgReC", "MiXCR", "pRESTO", "IgRepertoireConstructor"])
    plot_various_error_rate(igrec_dir + "/SIMULATED", what="prod", out="Fig_9c_prod_with_IGRC",
                            prod_criterion=True,
                            title="Sensitivity * precision (SIMULATED SIMPLE)",
                            kinds=["igrec", "mixcr2", "supernode", "ig_repertoire_constructor"], labels=["IgReC", "MiXCR", "pRESTO", "IgRepertoireConstructor"])



    plot_two_sums(igrec_dir + "/SIMULATED", out="Fig_10",
                  title="sensitivity + precision (SIMULATED)",
                  kind="igrec", label="IgReC")

    plot_two_sums(igrec_dir + "/SIMULATED", out="Fig_10_prod",
                  prod_criterion=True,
                  title="Sensitivity * precision (SIMULATED)",
                  kind="igrec", label="IgReC")

    plot_various_error_rate(igrec_dir + "/v_e_r_r", what="sensitivity", out="Fig_12_a",
                            title="Sensitivity (SYNTHETIC)",
                            kinds=["igrec", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])
    plot_various_error_rate(igrec_dir + "/v_e_r_r", what="precision", out="Fig_12_b",
                            title="Precision (SYNTHETIC)",
                            kinds=["igrec", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])

    plot_various_error_rate(igrec_dir + "/v_e_r_r", what="sensitivity", out="Fig_12_a_prod",
                            prod_criterion=True,
                            title="Sensitivity (SYNTHETIC)",
                            kinds=["igrec", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])
    plot_various_error_rate(igrec_dir + "/v_e_r_r", what="precision", out="Fig_12_b_prod",
                            prod_criterion=True,
                            title="Precision (SYNTHETIC)",
                            kinds=["igrec", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])


    sys.exit(0)
    plot_various_error_rate(igrec_dir + "/SIMULATED", out="sensitivity_plus_precision", what="sum",
                            title="SIMULATED, sensitivity + precision, simple",
                            kinds=["igrec", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])
    plot_various_error_rate(igrec_dir + "/SIMULATED", out="optimal_minsize", what="minsize",
                            title="SIMULATED, optimal min size, simple",
                            kinds=["igrec", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])


    plot_various_error_rate(igrec_dir + "/SIMULATED", what="sensitivity", out="sensitivity_woans", woans=True,
                            title="SIMULATED, sensitivity, complex",
                            kinds=["igrec_tau3", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])
    plot_various_error_rate(igrec_dir + "/SIMULATED", what="precision", out="precision_woans", woans=True,
                            title="SIMULATED, precision, complex",
                            kinds=["igrec_tau3", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])

    plot_various_error_rate(igrec_dir + "/v_e_r_r", out="sensitivity_plus_precision_flu", what="sum",
                            title="SYNTHETIC FLU, sensitivity + precision",
                            kinds=["igrec", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])

    # plot_various_error_rate(igrec_dir + "/SIMULATED",
    #                         kind1="igrec_tau3_f03", kind2="mixcr2",
    #                         what="precision", out="precision_tau3")
    # plot_various_error_rate(igrec_dir + "/SIMULATED",
    #                         kind1="igrec_tau3_f03", kind2="mixcr2",
    #                         what="sensitivity", out="sensitivity_tau3")
    rocs(igrec_dir + "/SIMULATED/errate_2.0000",
         tools=["igrec", "mixcr2", "supernode"],
         labels=["IgReC", "MiXCR", "pRESTO"],
         title="SIMULATED, error rate = 2, simple",
         out="sensitivity_precision_plot_er2")
    rocs(igrec_dir + "/SIMULATED/errate_1.0000",
         tools=["igrec", "mixcr2", "supernode"],
         labels=["IgReC", "MiXCR", "pRESTO"],
         title="SIMULATED, error rate = 1, simple",
         out="sensitivity_precision_plot_er1")
    rocs(igrec_dir + "/SIMULATED/errate_0.5000",
         tools=["igrec", "mixcr2", "supernode"],
         labels=["IgReC", "MiXCR", "pRESTO"],
         title="SIMULATED, error rate = 0.5, simple",
         out="sensitivity_precision_plot_er05")

    rocs(igrec_dir + "/SIMULATED/errate_2.0000_woans",
         tools=["igrec"],
         labels=["IgReC"],
         title="SIMULATED, error rate = 2, complex",
         out="sensitivity_precision_plot_er2_woanswer")
    rocs(igrec_dir + "/SIMULATED/errate_1.0000_woans",
         tools=["igrec"],
         labels=["IgReC"],
         title="SIMULATED, error rate = 1, complex",
         out="sensitivity_precision_plot_er1_woanswer")
    rocs(igrec_dir + "/SIMULATED/errate_0.5000_woans",
         tools=["igrec"],
         labels=["IgReC"],
         title="SIMULATED, error rate = 0.5, complex",
         out="sensitivity_precision_plot_er05_woanswer")

    rocs(igrec_dir + "/v_e_r_r/errate_2.0000",
         tools=["igrec", "mixcr2", "supernode"],
         labels=["IgReC", "MiXCR", "pRESTO"],
         title="SYNTHETIC, error rate = 2, simple",
         out="sensitivity_precision_plot_er2_real")
    rocs(igrec_dir + "/v_e_r_r/errate_1.0000",
         tools=["igrec", "mixcr2", "supernode"],
         labels=["IgReC", "MiXCR", "pRESTO"],
         title="SYNTHETIC, error rate = 1, simple",
         out="sensitivity_precision_plot_er1_real")
    rocs(igrec_dir + "/v_e_r_r/errate_0.5000",
         tools=["igrec", "mixcr2", "supernode"],
         labels=["IgReC", "MiXCR", "pRESTO"],
         title="SYNTHETIC, error rate = 0.5, simple",
         out="sensitivity_precision_plot_er05_real")
