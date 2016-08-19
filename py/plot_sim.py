#!/usr/bin/env python2

from plot_test import *

if __name__ == "__main__":
    rocs(igrec_dir + "/various_error_rate/errate_2.00",
         tools=["igrec_tau3", "mixcr", "supernode"],
         labels=["IgReC", "MiXCR", "pRESTO"],
         title="Error-rate = 2, with answer",
         out="sensitivity_precision_plot_er2")
    rocs(igrec_dir + "/various_error_rate/errate_1.00",
         tools=["igrec_tau3", "mixcr", "supernode"],
         labels=["IgReC", "MiXCR", "pRESTO"],
         title="Error-rate = 1, with answer",
         out="sensitivity_precision_plot_er1")
    rocs(igrec_dir + "/various_error_rate/errate_0.50",
         tools=["igrec_tau3", "mixcr", "supernode"],
         labels=["IgReC", "MiXCR", "pRESTO"],
         title="Error-rate = 0.5, with answer",
         out="sensitivity_precision_plot_er0.5")
    rocs(igrec_dir + "/various_error_rate/errate_2.00_woans",
         tools=["igrec_tau3"],
         labels=["IgReC"],
         title="Error-rate = 2, without answer",
         out="sensitivity_precision_plot_er2_woanswer")
    rocs(igrec_dir + "/various_error_rate/errate_1.00_woans",
         tools=["igrec_tau3"],
         labels=["IgReC"],
         title="Error-rate = 1, without answer",
         out="sensitivity_precision_plot_er1_woanswer")
    rocs(igrec_dir + "/various_error_rate/errate_0.50_woans",
         tools=["igrec_tau3"],
         labels=["IgReC"],
         title="Error-rate = 0.5, without answer",
         out="sensitivity_precision_plot_er0.5_woanswer")



    rocs(igrec_dir + "/v_e_r_r/errate_2.00",
         tools=["igrec_tau3", "supernode"],
         labels=["IgReC", "pRESTO"],
         title="Error-rate = 2, with answer",
         out="sensitivity_precision_plot_er2_real")
    rocs(igrec_dir + "/v_e_r_r/errate_1.00",
         tools=["igrec_tau3", "supernode"],
         labels=["IgReC", "pRESTO"],
         title="Error-rate = 1, with answer",
         out="sensitivity_precision_plot_er1_real")
    rocs(igrec_dir + "/v_e_r_r/errate_0.50",
         tools=["igrec_tau3", "supernode"],
         labels=["IgReC", "pRESTO"],
         title="Error-rate = 0.5, with answer",
         out="sensitivity_precision_plot_er0.5_real")


    plot_various_error_rate(igrec_dir + "/various_error_rate", what="sensitivity", out="sensitivity",
                            kinds=["igrec_tau3", "mixcr", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])
    plot_various_error_rate(igrec_dir + "/various_error_rate", what="precision", out="precision",
                            kinds=["igrec_tau3", "mixcr", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])
    plot_various_error_rate(igrec_dir + "/various_error_rate", what="sensitivity", out="sensitivity_woans", woans=True,
                            kinds=["igrec_tau3", "mixcr", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])
    plot_various_error_rate(igrec_dir + "/various_error_rate", what="precision", out="precision_woans", woans=True,
                            kinds=["igrec_tau3", "mixcr", "supernode"], labels=["IgReC", "MiXCR", "pRESTO"])

    plot_various_error_rate(igrec_dir + "/v_e_r_r", what="sensitivity", out="sensitivity_flu",
                            kinds=["igrec", "supernode"], labels=["IgReC", "pRESTO"])
    plot_various_error_rate(igrec_dir + "/v_e_r_r", what="precision", out="precision_flu",
                            kinds=["igrec", "supernode"], labels=["IgReC", "pRESTO"])
    # plot_various_error_rate(igrec_dir + "/various_error_rate",
    #                         kind1="igrec_tau3_f03", kind2="mixcr",
    #                         what="precision", out="precision_tau3")
    # plot_various_error_rate(igrec_dir + "/various_error_rate",
    #                         kind1="igrec_tau3_f03", kind2="mixcr",
    #                         what="sensitivity", out="sensitivity_tau3")
