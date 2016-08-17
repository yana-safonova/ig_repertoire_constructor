#!/usr/bin/env python2

from simulate import *
from aimquast_impl import plot_rocs
from aimquast_impl import plot_various_error_rate


# def two_rocs(dir, tool1, tool2, out, label1=None, label2=None, title=""):
#     if label1 is None:
#         label1 = tool1
#     if label2 is None:
#         label2 = tool2
#
#     plot_rocs([dir + "/" + tool1 + "/aimquast/aimquast.json",
#                dir + "/" + tool2 + "/aimquast/aimquast.json"],
#               labels=[label1, label2],
#               title=title,
#               format=("png",),
#               out=out)


def rocs(dir, tools, out, labels=None, title=""):
    if labels is None:
        labels = tools

    assert len(labels) == len(tools)

    plot_rocs([dir + "/" + tool + "/aimquast/aimquast.json" for tool in tools],
              labels=labels,
              title=title,
              format=("png",),
              out=out)


def two_rocs(dir, tool1, tool2, out, label1=None, label2=None, title=""):
    if label1 is None:
        label1 = tool1
    if label2 is None:
        label2 = tool2

    rocs(dir, [tool1, tool2], out, [label1, label2], title)


def plotplot(dir, out_dir, title):
    mkdir_p(out_dir)
    rocs(dir,
         tools=["igrec", "supernode"],
         labels=["IgReC", "pRESTO"],
         title="Real data: " + title,
         out=out_dir + "/sensitivity_precision_plot_all")
    rocs(dir,
         tools=["igrec_tau3", "supernode"],
         labels=["IgReC tau = 3", "pRESTO"],
         title="Real data: " + title,
         out=out_dir + "/sensitivity_precision_plot_all_tau3")
    rocs(dir,
         tools=["igrec_tau2", "supernode"],
         labels=["IgReC tau = 2", "pRESTO"],
         title="Real data: " + title,
         out=out_dir + "/sensitivity_precision_plot_all_tau2")
    rocs(dir,
         tools=["igrec_tau2", "supernode"],
         labels=["IgReC tau = 1", "pRESTO"],
         title="Real data: " + title,
         out=out_dir + "/sensitivity_precision_plot_all_tau1")
    #
    # two_rocs(dir,
    #          tool1="igrec", tool2="mixcr",
    #          label1="IgReC",
    #          label2="MiXCR",
    #          title="Real data: " + title,
    #          out=out_dir + "/sensitivity_precision_plot_igrec_mixcr")
    #
    # two_rocs(dir,
    #          tool1="igrec", tool2="supernode",
    #          label1="IgReC",
    #          label2="Supernode",
    #          title="Real data: " + title,
    #          out=out_dir + "/sensitivity_precision_plot_igrec_supernode")
    #
    # two_rocs(dir,
    #          tool1="igrec", tool2="igrec_tau1",
    #          label1="IgReC",
    #          label2="IgReC tau=1",
    #          title="Real data: " + title,
    #          out=out_dir + "/sensitivity_precision_plot_igrec_igrec_tau1")
    #
    # two_rocs(dir,
    #          tool1="igrec_tau1", tool2="supernode",
    #          label1="IgReC tau=1",
    #          label2="Supernode",
    #          title="Real data: " + title,
    #          out=out_dir + "/sensitivity_precision_plot_igrec_tau1_supernode")
    #
    # two_rocs(dir,
    #          tool1="mixcr", tool2="supernode",
    #          label1="MiXCR",
    #          label2="Supernode",
    #          title="Real data: " + title,
    #          out=out_dir + "/sensitivity_precision_plot_mixcr_supernode")
    #
    # two_rocs(dir,
    #          tool1="igrec_tau1_f09", tool2="supernode",
    #          label1="IgReC tau=1 f09",
    #          label2="Supernode",
    #          title="Real data: " + title,
    #          out=out_dir + "/sensitivity_precision_plot_igrec_tau1_f09_supernode")
    #
    # two_rocs(dir,
    #          tool1="igrec_tau2", tool2="supernode",
    #          label1="IgReC tau=2",
    #          label2="Supernode",
    #          title="Real data: " + title,
    #          out=out_dir + "/sensitivity_precision_plot_igrec_tau2_supernode")
    #
    # two_rocs(dir,
    #          tool1="igrec", tool2="igrec_tau2",
    #          label1="IgReC",
    #          label2="IgReC tau=2",
    #          title="Real data: " + title,
    #          out=out_dir + "/sensitivity_precision_plot_igrec_igrec_tau2")
    #
    # two_rocs(dir,
    #          tool1="igrec", tool2="igrec_tau3",
    #          label1="IgReC",
    #          label2="IgReC tau=3",
    #          title="Real data: " + title,
    #          out=out_dir + "/sensitivity_precision_plot_igrec_igrec_tau3")

if __name__ == "__main__":
    plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGH/filtering1/", "FLU_FV_21_IGH_1", title="FLU_FV_21_IGH_1")
    # plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGH/filtering2/", "FLU_FV_21_IGH_2", title="FLU_FV_21_IGH_2")
    # plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGH/filtering3/", "FLU_FV_21_IGH_3", title="FLU_FV_21_IGH_3")

    plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGL/filtering1/", "FLU_FV_21_IGL_1", title="FLU_FV_21_IGL_1")
    # plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGL/filtering2/", "FLU_FV_21_IGL_2", title="FLU_FV_21_IGL_2")
    # plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGL/filtering3/", "FLU_FV_21_IGL_3", title="FLU_FV_21_IGL_3")

    # plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGK/filtering1/", "FLU_FV_21_IGK_1", title="FLU_FV_21_IGK_1")
    # plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGK/filtering2/", "FLU_FV_21_IGK_2", title="FLU_FV_21_IGK_2")
    # plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGK/filtering3/", "FLU_FV_21_IGK_3", title="FLU_FV_21_IGK_3")

    plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGH/filtering1_jit0.5/", "FLU_FV_21_IGH_1_jit0.5", title="FLU_FV_21_IGH_1_jit0.5")
    plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGH/filtering1_jit1/", "FLU_FV_21_IGH_1_jit1", title="FLU_FV_21_IGH_1_jit1")
    plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGL/filtering1_jit0.5/", "FLU_FV_21_IGL_1_jit0.5", title="FLU_FV_21_IGL_1_jit0.5")
    plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGL/filtering1_jit1/", "FLU_FV_21_IGL_1_jit1", title="FLU_FV_21_IGL_1_jit1")

    sys.exit()

    plotplot(igrec_dir + "/src/extra/ig_quast_tool/AGE3/filtering3/", "AGE3_3", title="AGE3")
    plotplot(igrec_dir + "/src/extra/ig_quast_tool/AGE3/filtering2/", "AGE3_2", title="AGE3")
    plotplot(igrec_dir + "/src/extra/ig_quast_tool/AGE3/filtering1/", "AGE3_1", title="AGE3")
    plotplot(igrec_dir + "/src/extra/ref_bak_new/FLU_FV_27/filtering3/", "FLU_FV_27_3", title="FLU_FV_27")
    # plotplot(igrec_dir + "/src/extra/ref_bak_new/FLU_FV_27/filtering2/", "FLU_FV_27_2", title="FLU_FV_27")
    # plotplot(igrec_dir + "/src/extra/ref_bak_new/FLU_FV_27/filtering1/", "FLU_FV_27_1", title="FLU_FV_27")
    plotplot(igrec_dir + "/src/extra/ref_bak_new/FLU_FV_21/filtering3/", "FLU_FV_21_3", title="FLU_FV_21")
    # plotplot(igrec_dir + "/src/extra/ref_bak_new/FLU_FV_21/filtering2/", "FLU_FV_21_2", title="FLU_FV_21")
    # plotplot(igrec_dir + "/src/extra/ref_bak_new/FLU_FV_21/filtering1/", "FLU_FV_21_1", title="FLU_FV_21")

    plotplot(igrec_dir + "/src/extra/ref_bak_new/FLU_FV_22/filtering3/", "FLU_FV_22_3", title="FLU_FV_22")
    plotplot(igrec_dir + "/src/extra/ref_bak_new/FLU_FV_23/filtering3/", "FLU_FV_23_3", title="FLU_FV_23")
    sys.exit()

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
