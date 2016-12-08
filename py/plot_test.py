#!/usr/bin/env python2

from simulate import *
from aimquast_impl import *


def get_plot_various_error_rate_data(dir,
                                     kind="igrec",
                                     what="sensitivity",
                                     woans=False,
                                     size=5,
                                     add_aimquast_to_path=True):
    from glob import glob

    suff = "_woans" if woans else ""
    dirnames = glob(dir + "/errate_?.????" + suff)

    def extract_lambda(dirname):
        import re
        m = re.match(r".*/errate_([\d.]+)" + suff, dirname)
        if m:
            g = m.groups()
            return float(g[0].strip())
        else:
            return None

    if not woans:
        assert extract_lambda("fdsfsdfsd/errate_3.14") == 3.14

    def extract_data(dirname):
        jfile = "/aimquast/aimquast.json" if add_aimquast_to_path else "/aimquast.json"
        json = json_from_file(dirname + "/" + kind + jfile)
        return json

    lambdas = map(extract_lambda, dirnames)
    data = map(extract_data, dirnames)

    print lambdas
    lambdas, data = zip(*sorted([(l, d) for l, d in zip(lambdas, data) if l < 2.00001]))

    return lambdas, data


def get_plot_various_error_rate_data_serg(dir,
                                          kind="igrec",
                                          what="sensitivity",
                                          size=5,
                                          woans=False,
                                          add_aimquast_to_path=False):
    from glob import glob

    dirnames = glob(dir + "/pcr*")

    def extract_lambda(dirname):
        import re
        m = re.match(r".*pcr_([\d.]+)_super_100000_umi_15.*", dirname)
        if m:
            g = m.groups()
            return float(g[0].strip())
        else:
            return None

    if not woans:
        assert extract_lambda("pcr_0.003000_super_100000_umi_15") == 0.003

    def extract_data(dirname):
        jfile = "/quast_%s/aimquast.json" % kind
        json = json_from_file(dirname + "/" + jfile)
        return json

    nt_er_to_read_er = {0.0006: 0.534849, 0.0009 : 0.801606, 0.0012: 1.06461, 0.0015 : 1.33057, 0.0018: 1.60311, 0.0021 : 1.88114, 0.0024: 2.13592,
                        0.0025: 2.22422, 0.0027 : 2.41297, 0.0030: 2.66568, 0.0033 : 2.95098, 0.0036: 3.2174, 0.0039 : 3.48101, 0.0042: 3.75255,
                        0.0048: 4.29406, 0.0054: 4.81962, 0.006: 5.35227}
    lambdas = [nt_er_to_read_er[extract_lambda(dirname)] for dirname in dirnames]
    lambdas = [lam for lam in lambdas if lam <= 4]
    data = [extract_data(dirname) for dirname in dirnames if nt_er_to_read_er[extract_lambda(dirname)] in lambdas]

    lambdas, data = zip(*sorted([(l, d) for l, d in zip(lambdas, data) if l < 2000000000000.00001]))
    print lambdas

    return lambdas, data


def plot_various_error_rate_serg(dir,
                                 kinds,
                                 labels,
                                 sizes=None,
                                 out="var_error_rate",
                                 what="sensitivity", woans=False,
                                 title="",
                                 format=("png", "pdf", "svg"),
                                 legend=True,
                                 legend_loc=3,
                                 which=None):
    lambdas, _ = get_plot_various_error_rate_data_serg(dir, kind=kinds[0], woans=woans)
    datas = [get_plot_various_error_rate_data_serg(dir, kind=kind, woans=woans)[1] for kind in kinds]
    import matplotlib.pyplot as plt
    import seaborn as sns

    f, ax = initialize_plot()

    sns.set_style("darkgrid")
    colors = [tool2color(label) for label in labels]

    def opt_size(sensitivity, precision):
        return 1 + max(xrange(len(sensitivity)), key=lambda i: sensitivity[i] + precision[i])

    def get_what(dataset, what, cur_sizes):
        if what in ["sensitivity", "precision"]:
            return [data["reference_based"]["__data_" + what][size - 1] for data, size in zip(dataset, cur_sizes)]
        elif what == "sum":
            return [data["reference_based"]["__data_sensitivity"][size - 1] + data["reference_based"]["__data_precision"][size - 1] for data, size in zip(dataset, cur_sizes)]
        elif what == "minsize":
            return [size for data, size in zip(dataset, cur_sizes)]
        else:
            return None

    forplot = []
    for dataset in datas:
        if sizes is None:
            cur_sizes = [opt_size(data["reference_based"]["__data_sensitivity"], data["reference_based"]["__data_precision"]) for data in dataset]
        else:
            cur_sizes = [sizes] * len(dataset)
        cur_for_plot = get_what(dataset, what, cur_sizes)
        forplot.append(cur_for_plot)

    zipped = zip(forplot, colors, labels)
    if which is not None:
        zipped = [zipped[i] for i in which]
    for y, color, label in zipped:
        plt.plot(lambdas, y,
                 "b-", color=color, label=label)

    eps = 0.025
    if what in ["sensitivity", "precision"]:
        plt.ylim((0. - eps, 1. + eps))
    elif what == "sum":
        plt.ylim((1. - eps, 2. + eps))

    if what in ["sensitivity", "precision"]:
        plt.ylabel(what)
    elif what == "sum":
        plt.ylabel("sensitivity + precision")
    elif what == "minsize":
        plt.ylabel("optimal constructed min size")

    plt.xlabel("Error rate")

    if title:
        plt.title(title)

    if legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc=legend_loc)

    save_plot(out, format=format)




def plot_various_error_rate(dir,
                            kinds,
                            labels,
                            sizes=None,
                            out="var_error_rate",
                            what="sensitivity", woans=False,
                            title="",
                            format=("png", "pdf", "svg"),
                            legend=True,
                            which=None):
    lambdas, _ = get_plot_various_error_rate_data(dir, kind=kinds[0], woans=woans)
    datas = [get_plot_various_error_rate_data(dir, kind=kind, woans=woans)[1] for kind in kinds]
    import matplotlib.pyplot as plt
    import seaborn as sns

    f, ax = initialize_plot()

    sns.set_style("darkgrid")
    colors = [tool2color(label) for label in labels]

    def opt_size(sensitivity, precision):
        return 1 + max(xrange(len(sensitivity)), key=lambda i: sensitivity[i] + precision[i])

    def get_what(dataset, what, cur_sizes):
        if what in ["sensitivity", "precision"]:
            return [data["reference_based"]["__data_" + what][size - 1] for data, size in zip(dataset, cur_sizes)]
        elif what == "sum":
            return [data["reference_based"]["__data_sensitivity"][size - 1] + data["reference_based"]["__data_precision"][size - 1] for data, size in zip(dataset, cur_sizes)]
        elif what == "minsize":
            return [size for data, size in zip(dataset, cur_sizes)]
        else:
            return None

    forplot = []
    for dataset in datas:
        if sizes is None:
            cur_sizes = [opt_size(data["reference_based"]["__data_sensitivity"], data["reference_based"]["__data_precision"]) for data in dataset]
        else:
            cur_sizes = [sizes] * len(dataset)
        cur_for_plot = get_what(dataset, what, cur_sizes)
        forplot.append(cur_for_plot)

    zipped = zip(forplot, colors, labels)
    if which is not None:
        zipped = [zipped[i] for i in which]
    for y, color, label in zipped:
        plt.plot(lambdas, y,
                 "b-", color=color, label=label)

    eps = 0.025
    if what in ["sensitivity", "precision"]:
        plt.ylim((0. - eps, 1. + eps))
    elif what == "sum":
        plt.ylim((1. - eps, 2. + eps))

    if what in ["sensitivity", "precision"]:
        plt.ylabel(what)
    elif what == "sum":
        plt.ylabel("sensitivity + precision")
    elif what == "minsize":
        plt.ylabel("optimal constructed min size")

    plt.xlabel("Error rate")

    if title:
        plt.title(title)

    if legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc=3)

    save_plot(out, format=format)


def plot_two_sums(dir,
                  kind="igrec",
                  label="IgReC",
                  sizes=None,
                  out="var_error_rate",
                  title="",
                  legend=True,
                  format=("png", "pdf", "svg"),
                  which=None):
    lambdas, _ = get_plot_various_error_rate_data(dir, kind=kind, woans=False)
    data_wa = get_plot_various_error_rate_data(dir, kind=kind, woans=False)[1]
    data_woa = get_plot_various_error_rate_data(dir, kind=kind, woans=True)[1]
    import matplotlib.pyplot as plt
    import seaborn as sns

    f, ax = initialize_plot()

    sns.set_style("darkgrid")
    colors = ["cornflowerblue", "red", "orange", "black"]

    def opt_size(sensitivity, precision):
        return 1 + max(xrange(len(sensitivity)), key=lambda i: sensitivity[i] + precision[i])

    def get_what(dataset, what, cur_sizes):
        if what in ["sensitivity", "precision"]:
            return [data["reference_based"]["__data_" + what][size - 1] for data, size in zip(dataset, cur_sizes)]
        elif what == "sum":
            return [data["reference_based"]["__data_sensitivity"][size - 1] + data["reference_based"]["__data_precision"][size - 1] for data, size in zip(dataset, cur_sizes)]
        elif what == "minsize":
            return [size for data, size in zip(dataset, cur_sizes)]
        else:
            return None

    forplot = []
    for dataset in [data_wa, data_woa]:
        if sizes is None:
            cur_sizes = [opt_size(data["reference_based"]["__data_sensitivity"], data["reference_based"]["__data_precision"]) for data in dataset]
        else:
            cur_sizes = [sizes] * len(dataset)
        cur_for_plot = get_what(dataset, "sum", cur_sizes)
        forplot.append(cur_for_plot)

    labels = [label + ", simple", label + ", complex"]
    zipped = zip(forplot, colors, labels)
    if which is not None:
        zipped = [zipped[i] for i in which]
    for y, color, label in zipped:
        plt.plot(lambdas, y,
                 "b-", color=color, label=label)

    eps = 0.025
    plt.ylim((1. - eps, 2. + eps))

    plt.ylabel("sensitivity + precision")
    plt.xlabel("Error rate")

    if title:
        plt.title(title)

    if legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc=3)

    save_plot(out, format=format)


def tool2color(tool, secondary=False):
    primary_colors =   ["cornflowerblue", "seagreen", "orange", "black", "violet", "black"]
    secondary_colors = ["blue", "green", "darkorange", "dimgray", "orchid", "black"]
    colors = secondary_colors if secondary else primary_colors

    def tool2id(tool):
        tool = tool.lower()

        tools = ["igrec", "mixcr", "presto", "migec", "barigrec", "migec + mixcr"]
        from Levenshtein import distance
        dists = [distance(tool, t) for t in tools]
        return min(range(len(dists)), key=lambda i: dists[i])

    return colors[tool2id(tool)]


def plot_rocs(jsons, labels,
              max_size_threshold=75,
              out="two_rocs",
              title="",
              format=None,
              show_coords=True):
    import matplotlib.pyplot as plt
    import seaborn as sns

    jsons = [json for json in jsons if json is not None]

    how_many = len(jsons)
    jsons = map(json_from_file, jsons)

    sensitivities = [json["reference_based"]["__data_sensitivity"] for json in jsons]
    precisions = [json["reference_based"]["__data_precision"] for json in jsons]

    def opt_size(sensitivity, precision):
        return 1 + max(xrange(len(sensitivity)), key=lambda i: sensitivity[i] + precision[i])

    opt_sizes = [opt_size(sensitivity, precision) for sensitivity, precision in zip(sensitivities, precisions)]

    def cut_to_threshold(x):
        return x if len(x) <= max_size_threshold else x[:max_size_threshold]

    sensitivities = map(cut_to_threshold, sensitivities)
    precisions = map(cut_to_threshold, precisions)

    f, ax = initialize_plot()

    sns.set_style("darkgrid")

    # skip = 2
    skip = 0

    colors = [tool2color(label) for label in labels]
    point_colors = [tool2color(label, secondary=True) for label in labels]

    for precision, sensitivity, color, label in zip(precisions, sensitivities, colors, labels):
        plt.plot(precision[skip:], sensitivity[skip:], "b-", color=color, label=label)

    eps = 0.025
    plt.xlim((0 - eps, 1 + eps))
    plt.ylim((0 - eps, 1 + eps))

    plt.ylabel("sensitivity")
    plt.xlabel("precision")

    blue_points = []

    def annotation(size, x, y, color, text=None, xshift=0.01, yshift=-0.04,
                   add_text=False):
        _x = x[size - 1]
        _y = y[size - 1]

        bp = plt.plot(_x, _y, "bo", color=color, label="min size")
        blue_points.append(bp)

        if add_text:
            if text is None:
                text = str(size)
            plt.annotate(text, xy=(_x, _y),
                         color=color,
                         xytext=(_x + xshift, _y + yshift))

    for i in [1, 3, 5, 10, 50]:
        if i <= skip:
            continue
    # for i in []:
        annotation(i, precisions[0], sensitivities[0], color=point_colors[0])
        if len(jsons) > 1:
            annotation(i, precisions[1], sensitivities[1], color=point_colors[1], xshift=-0.04)
        if len(jsons) > 2:
            annotation(i, precisions[2], sensitivities[2], color=point_colors[2], yshift=0.04)
        if len(jsons) > 3:
            annotation(i, precisions[3], sensitivities[3], color=point_colors[3], xshift=-0.04, yshift=0.04)

    def red_point(x, y, size, label="NOLABEL", do_plot=True):
        if do_plot:
            plt.plot(x[size - 1], y[size - 1], "bo", color="red", label="Reference size threshold")
        if show_coords:
            print "%s, precision = %1.4f, sensitivity = %1.4f, minsize = %d" % (label, x[size - 1], y[size - 1], size)

    if show_coords:
        print "Dataset " + title
    for precision, sensitivity, color, label, opt_size in zip(precisions, sensitivities, colors, labels, opt_sizes):
        red_point(precision, sensitivity, opt_size, label)

    for precision, sensitivity, color, label, opt_size in zip(precisions, sensitivities, colors, labels, opt_sizes):
        red_point(precision, sensitivity, 5, label, do_plot=False)

    if title:
        plt.title(title)

    handles, labels = ax.get_legend_handles_labels()
    handles = handles[:how_many]
    labels = labels[:how_many]

    ax.legend(handles, labels, loc=3)

    save_plot(out, format=format)

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


def rocs(dir, tools, out, labels=None, title="", add_aimquast_to_path=True, **kwargs):
    if labels is None:
        labels = tools

    assert len(labels) == len(tools)

    jfile = "/aimquast/aimquast.json" if add_aimquast_to_path else "/aimquast.json"
    plot_rocs([dir + "/" + tool + jfile for tool in tools],
              labels=labels,
              title=title,
              format=("png", "pdf"),
              out=out,
              **kwargs)


def two_rocs(dir, tool1, tool2, out, label1=None, label2=None, title=""):
    if label1 is None:
        label1 = tool1
    if label2 is None:
        label2 = tool2

    rocs(dir, [tool1, tool2], out, [label1, label2], title)


def plotplot(dir, out_dir, title, **kwargs):
    mkdir_p(out_dir)
    rocs(dir,
         tools=["igrec", "mixcr2", "supernode"],
         labels=["IgReC", "MiXCR2", "pRESTO"],
         title=title,
         out=out_dir + "/sensitivity_precision_plot_all",
         **kwargs)
    rocs(dir,
         tools=["igrec", "mixcr", "supernode"],
         labels=["IgReC", "MiXCR", "pRESTO"],
         title=title,
         out=out_dir + "/sensitivity_precision_plot_all_old",
         **kwargs)
    rocs(dir,
         tools=["igrec_tau3", "mixcr", "supernode"],
         labels=["IgReC tau = 3", "MiXCR", "pRESTO"],
         title=title,
         out=out_dir + "/sensitivity_precision_plot_all_tau3",
         **kwargs)
    rocs(dir,
         tools=["igrec", "mixcr", "supernode", "igrec_vote"],
         labels=["IgReC", "MiXCR", "pRESTO", "IgReC split"],
         title=title,
         out=out_dir + "/sensitivity_precision_plot_all_split",
         **kwargs)
    # rocs(dir,
    #      tools=["igrec_tau2", "supernode"],
    #      labels=["IgReC tau = 2", "pRESTO"],
    #      title="Real data: " + title,
    #      out=out_dir + "/sensitivity_precision_plot_all_tau2")
    # rocs(dir,
    #      tools=["igrec_tau2", "supernode"],
    #      labels=["IgReC tau = 1", "pRESTO"],
    #      title="Real data: " + title,
    #      out=out_dir + "/sensitivity_precision_plot_all_tau1")


if __name__ == "__main__":
    plotplot(igrec_dir + "/src/extra/ig_quast_tool/AGE3/filtering3/", "AGE3_3", title="REAL", show_coords=True)
    # plotplot(igrec_dir + "/src/extra/ig_quast_tool/AGE3/filtering2/", "AGE3_2", title="AGE3")
    # plotplot(igrec_dir + "/src/extra/ig_quast_tool/AGE3/filtering1/", "AGE3_1", title="AGE3")

    plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGH/filtering3/", "FLU_FV_21_IGH_3", title="")
    sys.exit()
    plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGL/filtering3/", "FLU_FV_21_IGL_3", title="")
    plotplot(igrec_dir + "/src/extra/ig_quast_tool/AGE7/filtering3/", "AGE7_3", title="HEALTHY 2")


    plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGK/filtering3/", "FLU_FV_21_IGK_3", title="")

    plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGH/filtering3/", "FLU_FV_21_IGH_3", title="FLU_FV_21_IGH_3")
    plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_27_IGH/filtering3/", "FLU_FV_27_IGH_3", title="FLU_FV_27_IGH_3")

    plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGL/filtering3/", "FLU_FV_21_IGL_3", title="FLU_FV_21_IGL_3")
    plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_27_IGL/filtering3/", "FLU_FV_27_IGL_3", title="FLU_FV_27_IGL_3")

    plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_21_IGK/filtering3/", "FLU_FV_21_IGK_3", title="FLU_FV_21_IGK_3")
    plotplot(igrec_dir + "/src/extra/ig_quast_tool/FLU_FV_27_IGK/filtering3/", "FLU_FV_27_IGK_3", title="FLU_FV_27_IGK_3")

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
