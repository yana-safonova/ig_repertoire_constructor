#!/usr/bin/env python2

from simulate import *
from aimquast_impl import plot_two_rocs
from aimquast_impl import plot_various_error_rate


if __name__ == "__main__":
    plot_two_rocs(igrec_dir + "/src/extra/ig_quast_tool/AGE3/filtering3/igrec/aimquast/aimquast.json",
                  igrec_dir + "/src/extra/ig_quast_tool/AGE3/filtering3/mixcr/aimquast/aimquast.json",
                  label1="IgReC",
                  label2="MiXCR",
                  title="Real data: AGE3",
                  out="sensitivity_precision_plot_age3")

    plot_two_rocs(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_27/filtering3/igrec/aimquast/aimquast.json",
                  igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_27/filtering3/mixcr/aimquast/aimquast.json",
                  label1="IgReC",
                  label2="MiXCR",
                  title="Real data: FLU_FV_27",
                  out="sensitivity_precision_plot_flu_fv_27")

    plot_two_rocs(igrec_dir + "/various_error_rate/errate_2.00/igrec/aimquast/aimquast.json",
                  label1="IgReC",
                  label2="MiXCR",
                  title="Error-rate = 2, with answer",
                  out="sensitivity_precision_plot_er2")

    plot_two_rocs(igrec_dir + "/various_error_rate/errate_1.00/igrec_tau3_f03/aimquast/aimquast.json",
                  igrec_dir + "/various_error_rate/errate_1.00/mixcr/aimquast/aimquast.json",
                  label1="IgReC",
                  label2="MiXCR",
                  title="Error-rate = 1, with answer",
                  out="sensitivity_precision_plot_er1")

    plot_two_rocs(igrec_dir + "/various_error_rate/errate_2.00_woans/igrec_tau3_f03/aimquast/aimquast.json",
                  igrec_dir + "/various_error_rate/errate_2.00_woans/mixcr/aimquast/aimquast.json",
                  label1="IgReC",
                  label2="MiXCR",
                  plot_second=False,
                  title="Error-rate = 2, without answer",
                  out="sensitivity_precision_plot_er2_woanswer")

    plot_two_rocs(igrec_dir + "/various_error_rate/errate_1.00_woans/igrec_tau3_f03/aimquast/aimquast.json",
                  igrec_dir + "/various_error_rate/errate_1.00_woans/mixcr/aimquast/aimquast.json",
                  label1="IgReC",
                  label2="MiXCR",
                  plot_second=False,
                  title="Error-rate = 1, without answer",
                  out="sensitivity_precision_plot_er1_woanswer")

    plot_two_rocs(igrec_dir + "/various_error_rate/errate_0.50/igrec/aimquast/aimquast.json",
                  igrec_dir + "/various_error_rate/errate_0.50/mixcr/aimquast/aimquast.json",
                  label1="IgReC",
                  label2="MiXCR",
                  title="Error-rate = 0.5",
                  out="sensitivity_precision_plot_er05")

    plot_various_error_rate(igrec_dir + "/various_error_rate", what="sensitivity", out="sensitivity", kind1="igrec_tau3_f03", label1="IgReC", label2="MiXCR")
    plot_various_error_rate(igrec_dir + "/various_error_rate", what="precision", out="precision")
    plot_various_error_rate(igrec_dir + "/various_error_rate", what="sensitivity", out="sensitivity_woans", woans=True, kind1="igrec_tau3_f03", label1="IgReC", label2="MiXCR")
    plot_various_error_rate(igrec_dir + "/various_error_rate", what="precision", out="precision_woans", woans=True)

    plot_two_rocs(igrec_dir + "/various_error_rate/errate_2.00/igrec_tau3_f03/aimquast/aimquast.json",
                  igrec_dir + "/various_error_rate/errate_2.00/mixcr/aimquast/aimquast.json",
                  out="ROC_tau3")
    plot_various_error_rate(igrec_dir + "/various_error_rate",
                            kind1="igrec_tau3_f03", kind2="mixcr",
                            what="precision", out="precision_tau3")
    plot_various_error_rate(igrec_dir + "/various_error_rate",
                            kind1="igrec_tau3_f03", kind2="mixcr",
                            what="sensitivity", out="sensitivity_tau3")
