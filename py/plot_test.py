#!/usr/bin/env python2

from simulate import *
from aimquast_impl import plot_two_rocs
from aimquast_impl import plot_various_error_rate


def two_rocs(dir, tool1, tool2, out, label1=None, label2=None, title=""):
    if label1 is None:
        label1 = tool1
    if label2 is None:
        label2 = tool2

    plot_two_rocs(dir + "/" + tool1 + "/aimquast/aimquast.json",
                  dir + "/" + tool2 + "/aimquast/aimquast.json",
                  label1=label1, label2=label2,
                  title=title,
                  out=out)


if __name__ == "__main__":
    two_rocs(igrec_dir + "/src/extra/ig_quast_tool/AGE3/filtering3/",
             tool1="igrec", tool2="mixcr",
             label1="IgReC",
             label2="MiXCR",
             title="Real data: AGE3",
             out="sensitivity_precision_plot_age3")

    two_rocs(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21/filtering3/",
             tool1="igrec", tool2="mixcr",
             label1="IgReC",
             label2="MiXCR",
             title="Real data: FLU_FV_21",
             out="sensitivity_precision_plot_flu_fv_21")

    two_rocs(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21/filtering3/",
             tool1="igrec", tool2="supernode",
             label1="IgReC",
             label2="Supernode",
             title="Real data: FLU_FV_21",
             out="sensitivity_precision_plot_flu_fv_21")

    two_rocs(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_22/filtering3/",
             tool1="igrec", tool2="mixcr",
             label1="IgReC",
             label2="MiXCR",
             title="Real data: FLU_FV_22",
             out="sensitivity_precision_plot_flu_fv_22")

    two_rocs(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_27/filtering3/",
             tool1="igrec", tool2="mixcr",
             label1="IgReC",
             label2="MiXCR",
             title="Real data: FLU_FV_27",
             out="sensitivity_precision_plot_flu_fv_27")

    two_rocs(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_27/filtering3/",
             tool1="igrec", tool2="supernode",
             label1="IgReC",
             label2="Supernode",
             title="Real data: FLU_FV_27",
             out="sensitivity_precision_plot_flu_fv_27")

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

    plot_various_error_rate(igrec_dir + "/various_error_rate", what="sensitivity", out="sensitivity", kind1="igrec_tau3_f03", label1="IgReC", label2="MiXCR")
    plot_various_error_rate(igrec_dir + "/various_error_rate", what="precision", out="precision")
    plot_various_error_rate(igrec_dir + "/various_error_rate", what="sensitivity", out="sensitivity_woans", woans=True, kind1="igrec_tau3_f03", label1="IgReC", label2="MiXCR")
    plot_various_error_rate(igrec_dir + "/various_error_rate", what="precision", out="precision_woans", woans=True)

    two_rocs(igrec_dir + "/various_error_rate/errate_2.00",
             tool1="igrec", tool2="mixcr",
             out="ROC_tau3")
    plot_various_error_rate(igrec_dir + "/various_error_rate",
                            kind1="igrec_tau3_f03", kind2="mixcr",
                            what="precision", out="precision_tau3")
    plot_various_error_rate(igrec_dir + "/various_error_rate",
                            kind1="igrec_tau3_f03", kind2="mixcr",
                            what="sensitivity", out="sensitivity_tau3")
