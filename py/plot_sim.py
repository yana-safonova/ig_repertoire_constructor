#!/usr/bin/env python2

from plot_test import *

if __name__ == "__main__":
    plot_various_error_rate(igrec_dir + "/SIMULATED", what="sensitivity", out="TEST_Fig_9a",
                            title="Sensitivity (SIMULATEDx10 SIMPLE)",
                            kinds=["igrec", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"],
                            multiple=True)
    plot_various_error_rate(igrec_dir + "/SIMULATED", what="precision", out="TEST_Fig_9b",
                            title="Precision (SIMULATEDx10 SIMPLE)",
                            kinds=["igrec", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"],
                            multiple=True)
    plot_various_error_rate(igrec_dir + "/SIMULATED", what="sum", out="TEST_Fig_9c",
                            title="Sensitivity + precision (SIMULATEDx10 SIMPLE)",
                            kinds=["igrec", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"],
                            multiple=True)
    plot_various_error_rate(igrec_dir + "/SIMULATED", what="minsize", out="TEST_Fig_9d",
                            title="Optimal minsize (SIMULATEDx10 SIMPLE)",
                            kinds=["igrec", "mixcr2", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"],
                            which=[1, 2],  # remove IgReC's line
                            multiple=True)

    sys.exit()
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
