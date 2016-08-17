#!/usr/bin/env python2

from plot_test import *

if __name__ == "__main__":
    two_rocs(igrec_dir + "/various_error_rate/errate_2.00",
             tool1="igrec", tool2="mixcr",
             label1="IgReC",
             label2="MiXCR",
             title="Error-rate = 2, with answer",
             out="sensitivity_precision_plot_er2")

    two_rocs(igrec_dir + "/various_error_rate/errate_1.00",
             tool1="igrec", tool2="mixcr",
             label1="IgReC",
             label2="MiXCR",
             title="Error-rate = 1, with answer",
             out="sensitivity_precision_plot_er1")

    two_rocs(igrec_dir + "/various_error_rate/errate_2.00_woans",
             tool1="igrec", tool2="mixcr",
             label1="IgReC",
             label2="MiXCR",
             plot_second=False,
             title="Error-rate = 2, without answer",
             out="sensitivity_precision_plot_er2_woanswer")

    two_rocs(igrec_dir + "/various_error_rate/errate_1.00_woans",
             tool1="igrec", tool2="mixcr",
             label1="IgReC",
             label2="MiXCR",
             plot_second=False,
             title="Error-rate = 1, without answer",
             out="sensitivity_precision_plot_er1_woanswer")

    two_rocs(igrec_dir + "/various_error_rate/errate_0.50",
             tool1="igrec", tool2="mixcr",
             label1="IgReC",
             label2="MiXCR",
             title="Error-rate = 0.5",
             out="sensitivity_precision_plot_er05")

    plot_various_error_rate(igrec_dir + "/various_error_rate", what="sensitivity", out="sensitivity",
                            kinds=["igrec_tau3_f03", "mixcr", "pRESTO"], labels=["IgReC", "MiXCR", "supernode"])
    plot_various_error_rate(igrec_dir + "/various_error_rate", what="precision", out="precision",
                            kinds=["igrec_tau3_f03", "mixcr", "pRESTO"], labels=["IgReC", "MiXCR", "supernode"])
    plot_various_error_rate(igrec_dir + "/various_error_rate", what="sensitivity", out="sensitivity_woans", woans=True,
                            kinds=["igrec_tau3_f03", "mixcr", "pRESTO"], labels=["IgReC", "MiXCR", "supernode"])
    plot_various_error_rate(igrec_dir + "/various_error_rate", what="precision", out="precision_woans", woans=True,
                            kinds=["igrec_tau3_f03", "mixcr", "pRESTO"], labels=["IgReC", "MiXCR", "supernode"])


    two_rocs(igrec_dir + "/various_error_rate/errate_2.00",
             tool1="igrec", tool2="mixcr",
             out="ROC_tau3")
    plot_various_error_rate(igrec_dir + "/various_error_rate",
                            kind1="igrec_tau3_f03", kind2="mixcr",
                            what="precision", out="precision_tau3")
    plot_various_error_rate(igrec_dir + "/various_error_rate",
                            kind1="igrec_tau3_f03", kind2="mixcr",
                            what="sensitivity", out="sensitivity_tau3")
