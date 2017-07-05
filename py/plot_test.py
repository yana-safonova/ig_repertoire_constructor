#!/usr/bin/env python2

from simulate import *
from igquast_impl import *


def get_plot_various_error_rate_data(dir,
                                     kind="igrec",
                                     what="sensitivity",
                                     woans=False,
                                     size=5,
                                     add_aimquast_to_path=True,
                                     multiple=False):
    from glob import glob

    suff = "_woans" if woans else ""
    mul_suff = "_seed_*" if multiple else ""
    dirnames = glob(dir + "/errate_?.????" + suff + mul_suff)

    mul_suff = r"_seed_\d+" if multiple else ""

    def extract_lambda(dirname):
        import re
        m = re.match(r".*/errate_([\d.]+)" + suff + mul_suff, dirname)
        if m:
            g = m.groups()
            return float(g[0].strip())
        else:
            return None

    if not woans and not multiple:
        assert extract_lambda("fdsfsdfsd/errate_3.14") == 3.14

    def extract_data(dirname):
        jfile = "/aimquast/aimquast.json" if add_aimquast_to_path else "/aimquast.json"
        json = json_from_file(dirname + "/" + kind + jfile)
        return json

    def unique(x):
        return sorted(list(set(x)))

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
                                 which=None,
                                 prod_criterion=False):
    lambdas, _ = get_plot_various_error_rate_data_serg(dir, kind=kinds[0], woans=woans)
    datas = [get_plot_various_error_rate_data_serg(dir, kind=kind, woans=woans)[1] for kind in kinds]
    import matplotlib.pyplot as plt
    import seaborn as sns

    f, ax = initialize_plot()

    sns.set_style("darkgrid")
    colors = [tool2color(label) for label in labels]

    if prod_criterion:
        def opt_size(sensitivity, precision):
            return 1 + max(xrange(len(sensitivity)), key=lambda i: sensitivity[i] * precision[i])
    else:
        def opt_size(sensitivity, precision):
            return 1 + max(xrange(len(sensitivity)), key=lambda i: sensitivity[i] + precision[i])

    def get_what(dataset, what, cur_sizes):
        if what in ["sensitivity", "precision"]:
            return [data["reference_based"]["__data_" + what][size - 1] for data, size in zip(dataset, cur_sizes)]
        elif what == "sum":
            return [data["reference_based"]["__data_sensitivity"][size - 1] + data["reference_based"]["__data_precision"][size - 1] for data, size in zip(dataset, cur_sizes)]
        elif what == "prod":
            return [data["reference_based"]["__data_sensitivity"][size - 1] * data["reference_based"]["__data_precision"][size - 1] for data, size in zip(dataset, cur_sizes)]
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

    print title, what
    for y, color, label in zipped:
        plt.plot(lambdas, y,
                 "b-", color=color, label=label)
        print label, y

    eps = 0.025
    if what in ["sensitivity", "precision"]:
        plt.ylim((0. - eps, 1. + eps))
    elif what == "sum":
        plt.ylim((1. - eps, 2. + eps))
    elif what == "prod":
        plt.ylim((0. - eps, 1. + eps))

    if what in ["sensitivity", "precision"]:
        plt.ylabel(what)
    elif what == "sum":
        plt.ylabel("sensitivity + precision")
    elif what == "prod":
        plt.ylabel("sensitivity * precision")
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
                            which=None,
                            prod_criterion=False,
                            multiple=False):
    lambdas, _ = get_plot_various_error_rate_data(dir, kind=kinds[0], woans=woans, multiple=multiple)
    datas = [get_plot_various_error_rate_data(dir, kind=kind, woans=woans, multiple=multiple)[1] for kind in kinds]
    import matplotlib.pyplot as plt
    import seaborn as sns

    f, ax = initialize_plot()

    sns.set_style("darkgrid")
    colors = [tool2color(label) for label in labels]

    if prod_criterion:
        def opt_size(sensitivity, precision):
            return 1 + max(xrange(len(sensitivity)), key=lambda i: sensitivity[i] * precision[i])
    else:
        def opt_size(sensitivity, precision):
            return 1 + max(xrange(len(sensitivity)), key=lambda i: sensitivity[i] + precision[i])

    def get_what(dataset, what, cur_sizes):
        if what in ["sensitivity", "precision"]:
            return [data["reference_based"]["__data_" + what][size - 1] for data, size in zip(dataset, cur_sizes)]
        elif what == "sum":
            return [data["reference_based"]["__data_sensitivity"][size - 1] + data["reference_based"]["__data_precision"][size - 1] for data, size in zip(dataset, cur_sizes)]
        elif what == "prod":
            return [data["reference_based"]["__data_sensitivity"][size - 1] * data["reference_based"]["__data_precision"][size - 1] for data, size in zip(dataset, cur_sizes)]
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

    print title, what
    for y, color, label in zipped:
        min_lambda = min(lambdas)
        nseeds = sum(1 for l in lambdas if l == min_lambda)
        nlambdas = len(lambdas) / nseeds
        if what == "minsize":
            means = []
            lms = []
            for i in range(nlambdas):
                _y = y[i*nseeds:(i+1)*nseeds]
                means.append(float(sum(_y)) / len(_y))
                lms.append(lambdas[i*nseeds])
                plt.plot([lambdas[i*nseeds]] * 2, [min(_y), max(_y)],
                         "--bo", color=color, label=label)
            plt.plot(lms, means,
                     "b-", color=color, label=label)
        else:
            for seed in range(nseeds):
                plt.plot(lambdas[seed::nseeds], y[seed::nseeds],
                         "b-", color=color, label=label)
        print label, y

    eps = 0.025
    if what in ["sensitivity", "precision"]:
        plt.ylim((0. - eps, 1. + eps))
    elif what == "sum":
        plt.ylim((1. - eps, 2. + eps))
    elif what == "prod":
        plt.ylim((0. - eps, 1. + eps))

    plt.xlim(((0. - eps) * max(lambdas), (1. + eps) * max(lambdas)))

    if what in ["sensitivity", "precision"]:
        plt.ylabel(what)
    elif what == "sum":
        plt.ylabel("sensitivity + precision")
    elif what == "prod":
        plt.ylabel("sensitivity * precision")
    elif what == "minsize":
        plt.ylabel("optimal constructed min size")

    plt.xlabel("Error rate")

    if title:
        plt.title(title)

    if legend:
        handles, labels = ax.get_legend_handles_labels()
        nlines = len(handles) / len(zipped)
        ax.legend(handles[nlines-1::nlines], labels[nlines-1::nlines], loc=3)

    save_plot(out, format=format)


def plot_two_sums(dir,
                  kind="igrec",
                  label="IgReC",
                  sizes=None,
                  out="var_error_rate",
                  title="",
                  legend=True,
                  format=("png", "pdf", "svg"),
                  which=None,
                  prod_criterion=False):
    lambdas, _ = get_plot_various_error_rate_data(dir, kind=kind, woans=False)
    data_wa = get_plot_various_error_rate_data(dir, kind=kind, woans=False)[1]
    data_woa = get_plot_various_error_rate_data(dir, kind=kind, woans=True)[1]
    import matplotlib.pyplot as plt
    import seaborn as sns

    f, ax = initialize_plot()

    sns.set_style("darkgrid")
    colors = ["cornflowerblue", "red", "orange", "black"]

    if prod_criterion:
        def opt_size(sensitivity, precision):
            return 1 + max(xrange(len(sensitivity)), key=lambda i: sensitivity[i] * precision[i])
        what = "prod"
    else:
        def opt_size(sensitivity, precision):
            return 1 + max(xrange(len(sensitivity)), key=lambda i: sensitivity[i] + precision[i])
        what = "sum"

    def get_what(dataset, what, cur_sizes):
        if what in ["sensitivity", "precision"]:
            return [data["reference_based"]["__data_" + what][size - 1] for data, size in zip(dataset, cur_sizes)]
        elif what == "sum":
            return [data["reference_based"]["__data_sensitivity"][size - 1] + data["reference_based"]["__data_precision"][size - 1] for data, size in zip(dataset, cur_sizes)]
        elif what == "prod":
            return [data["reference_based"]["__data_sensitivity"][size - 1] * data["reference_based"]["__data_precision"][size - 1] for data, size in zip(dataset, cur_sizes)]
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

        cur_for_plot = get_what(dataset, what, cur_sizes)
        forplot.append(cur_for_plot)

    labels = [label + ", SIMULATED SIMPLE", label + ", SIMULATED COMPLEX"]
    zipped = zip(forplot, colors, labels)
    if which is not None:
        zipped = [zipped[i] for i in which]
    for y, color, label in zipped:
        plt.plot(lambdas, y,
                 "b-", color=color, label=label)

    eps = 0.025
    if what == "sum":
        plt.ylim((1. - eps, 2. + eps))
        plt.ylabel("sensitivity + precision")
    elif what == "prod":
        plt.ylim((0. - eps, 1. + eps))
        plt.ylabel("sensitivity * precision")

    plt.xlabel("Error rate")

    if title:
        plt.title(title)

    if legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc=3)

    save_plot(out, format=format)


def tool2color(tool, secondary=False):
    primary_colors =   ["cornflowerblue", "seagreen", "orange", "black", "violet", "black", "lightpink", "gold"]
    secondary_colors = ["blue", "green", "darkorange", "dimgray", "orchid", "black", "hotpink", "yellow"]
    colors = secondary_colors if secondary else primary_colors

    def tool2id(tool):
        tool = tool.lower()

        tools = ["igrec", "mixcr", "presto", "migec", "barigrec", "migec + mixcr", "ig_repertoire_constructor", "igrec_tau3"]
        from Levenshtein import distance
        dists = [distance(tool, t) for t in tools]
        return min(range(len(dists)), key=lambda i: dists[i])

    return colors[tool2id(tool)]


def plot_rocs_seed(jsons, labels,
                   max_size_threshold=75,
                   out="two_rocs",
                   title="",
                   format=None,
                   show_coords=True,
                   prod_criterion=False):
    import matplotlib.pyplot as plt
    import seaborn as sns

    jsons = [json for json in jsons if json is not None]

    f, ax = initialize_plot()

    sns.set_style("darkgrid")

    # skip = 2
    skip = 0

    how_many = len(jsons)
    pat_jsons = jsons
    for seed in range(10):
        jsons = [json.replace("%%SEED%%", str(seed)) for json in pat_jsons]
        try:
            jsons = map(json_from_file, jsons)
        except:
            continue

        sensitivities = [json["reference_based"]["__data_sensitivity"] for json in jsons]
        precisions = [json["reference_based"]["__data_precision"] for json in jsons]

        if prod_criterion:
            def opt_size(sensitivity, precision):
                return 1 + max(xrange(len(sensitivity)), key=lambda i: sensitivity[i] * precision[i])
        else:
            def opt_size(sensitivity, precision):
                return 1 + max(xrange(len(sensitivity)), key=lambda i: sensitivity[i] + precision[i])

        opt_sizes = [opt_size(sensitivity, precision) for sensitivity, precision in zip(sensitivities, precisions)]

        def cut_to_threshold(x):
            return x if len(x) <= max_size_threshold else x[:max_size_threshold]

        sensitivities = map(cut_to_threshold, sensitivities)
        precisions = map(cut_to_threshold, precisions)

        colors = [tool2color(label) for label in labels]
        point_colors = [tool2color(label, secondary=True) for label in labels]

        for precision, sensitivity, color, label in zip(precisions, sensitivities, colors, labels):
            plt.plot(precision[skip:], sensitivity[skip:], "b-", color=color, label=label)

    eps = 0.025
    plt.xlim((0 - eps, 1 + eps))
    plt.ylim((0 - eps, 1 + eps))

    plt.ylabel("sensitivity")
    plt.xlabel("precision")

    if title:
        plt.title(title)

    handles, labels = ax.get_legend_handles_labels()
    handles = handles[:how_many]
    labels = labels[:how_many]

    ax.legend(handles, labels, loc=3)

    save_plot(out, format=format)


def remove_trailing_zeros(sensitivity, precision):
    assert len(sensitivity) == len(precision)
    # TODO add this to igQIUAST
    while (sensitivity and precision and sensitivity[-1] == precision[-1] == 0):
        # print "Trailing zero cropped"
        sensitivity = sensitivity[:-1]
        precision = precision[:-1]
    # print "All zeros cropped"
    assert len(sensitivity) == len(precision)
    return sensitivity, precision


def plot_rocs(jsons, labels,
              max_size_threshold=75,
              out="two_rocs",
              title="",
              format=None,
              show_coords=True,
              prod_criterion=False):
    import matplotlib.pyplot as plt
    import seaborn as sns

    jsons = [json for json in jsons if json is not None]

    how_many = len(jsons)
    jsons = map(json_from_file, jsons)

    sensitivities = [json["reference_based"]["__data_sensitivity"] for json in jsons]
    precisions = [json["reference_based"]["__data_precision"] for json in jsons]

    if prod_criterion:
        def opt_size(sensitivity, precision):
            return 1 + max(xrange(len(sensitivity)), key=lambda i: sensitivity[i] * precision[i])
    else:
        def opt_size(sensitivity, precision):
            return 1 + max(xrange(len(sensitivity)), key=lambda i: sensitivity[i] + precision[i])

    opt_sizes = [opt_size(sensitivity, precision) for sensitivity, precision in zip(sensitivities, precisions)]

    def cut_to_threshold(x):
        return x if len(x) <= max_size_threshold else x[:max_size_threshold]

    f, ax = initialize_plot()

    sns.set_style("darkgrid")

    # skip = 2
    skip = 0

    sensitivities = map(cut_to_threshold, sensitivities)
    precisions = map(cut_to_threshold, precisions)

    sensitivities, precisions = zip(*map(lambda x: remove_trailing_zeros(*x), zip(sensitivities, precisions)))

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
        if (len(x) < size or len(y) < size):
            return
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
    if "%%SEED%%" in dir:
        f = plot_rocs_seed
    else:
        f = plot_rocs
    f([dir + "/" + tool + jfile for tool in tools],
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
         labels=["IgReC", "MiXCR", "pRESTO"],
         title=title,
         out=out_dir + "/sensitivity_precision_plot",
         **kwargs)

    return
    rocs(dir,
         tools=["igrec", "mixcr2full", "supernode"],
         labels=["IgReC", "MiXCR_full", "pRESTO"],
         title=title,
         out=out_dir + "/sensitivity_precision_plot_full_mixcr",
         **kwargs)
    try:
        rocs(dir,
             tools=["igrec", "ig_repertoire_constructor", "igrec_tau3"],
             labels=["IgReC", "IgRepertoireConstructor", "IgReC tau=3"],
             title=title,
             out=out_dir + "/sensitivity_precision_plot_old_vs_new",
             **kwargs)
    except:
        pass

    try:
        rocs(dir,
             tools=["igrec", "mixcr2", "supernode", "ig_repertoire_constructor"],
             labels=["IgReC", "MiXCR", "pRESTO", "IgRepertoireConstructor"],
             title=title,
             out=out_dir + "/sensitivity_precision_plot_all",
             **kwargs)
    except:
        pass
    # rocs(dir,
    #      tools=["igrec", "mixcr", "supernode"],
    #      labels=["IgReC", "MiXCR", "pRESTO"],
    #      title=title,
    #      out=out_dir + "/sensitivity_precision_plot_all_old",
    #      **kwargs)
    # rocs(dir,
    #      tools=["igrec_tau3", "mixcr", "supernode"],
    #      labels=["IgReC tau = 3", "MiXCR", "pRESTO"],
    #      title=title,
    #      out=out_dir + "/sensitivity_precision_plot_all_tau3",
    #      **kwargs)
    # rocs(dir,
    #      tools=["igrec", "mixcr", "supernode", "igrec_vote"],
    #      labels=["IgReC", "MiXCR", "pRESTO", "IgReC split"],
    #      title=title,
    #      out=out_dir + "/sensitivity_precision_plot_all_split",
    #      **kwargs)
    # rocs(dir,
    #      tools=["igrec", "mixcr", "supernode"],
    #      labels=["IgReC", "MiXCR1", "pRESTO"],
    #      title=title,
    #      out=out_dir + "/sensitivity_precision_plot_all_old",
    #      **kwargs)
    # rocs(dir,
    #      tools=["igrec_tau3", "mixcr", "supernode"],
    #      labels=["IgReC tau = 3", "MiXCR1", "pRESTO"],
    #      title=title,
    #      out=out_dir + "/sensitivity_precision_plot_all_tau3",
    #      **kwargs)
    # rocs(dir,
    #      tools=["igrec", "mixcr", "supernode", "igrec_vote"],
    #      labels=["IgReC", "MiXCR", "pRESTO", "IgReC split"],
    #      title=title,
    #      out=out_dir + "/sensitivity_precision_plot_all_split",
    #      **kwargs)
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
    plotplot(igrec_dir + "py/test_on_pd/SIMTCR_0.5/", "SIMTCR_0.5_figs", title="SIM TCR SIMPLE dataset, 0.5 errors per read", show_coords=True)
    plotplot(igrec_dir + "py/test_on_pd/SIMTCR_1/", "SIMTCR_1_figs", title="SIM TCR SIMPLE dataset, 1 error per read", show_coords=True)
    plotplot(igrec_dir + "py/test_on_pd/SIMTCR_2/", "SIMTCR_2_figs", title="SIM TCR SIMPLE dataset, 2 errors per read", show_coords=True)
    plotplot(igrec_dir + "SIMULATED/errate_0.5000_seed_%%SEED%%", "MULT_SIMULATED_0.5_figs", title="SIMULATEDx10 SIMPLE dataset, 0.5 error per read", show_coords=True)
    plotplot(igrec_dir + "SIMULATED/errate_1.0000_seed_%%SEED%%", "MULT_SIMULATED_1_figs", title="SIMULATEDx10 SIMPLE dataset, 1 error per read", show_coords=True)
    plotplot(igrec_dir + "SIMULATED/errate_2.0000_seed_%%SEED%%", "MULT_SIMULATED_2_figs", title="SIMULATEDx10 SIMPLE dataset, 2 error per read", show_coords=True)

    plotplot(igrec_dir + "py/test_on_pd/REAL/", "REAL_figs", title="Sensitivity-precision plot (REAL dataset)", show_coords=True)
    plotplot(igrec_dir + "py/test_on_pd/REAL_CHU/", "REAL_CHU_figs", title="Sensitivity-precision plot (REAL MiGEC dataset)", show_coords=True)
    plotplot(igrec_dir + "py/test_on_pd/SIMULATED_1/", "SIMULATED_1_figs", title="SIMULATED SIMPLE dataset, 1 error per read", show_coords=True)
    plotplot(igrec_dir + "py/test_on_pd/SIMULATED_2/", "SIMULATED_2_figs", title="SIMULATED SIMPLE dataset, 2 errors per read", show_coords=True)
    plotplot(igrec_dir + "py/test_on_pd/SIMULATED_0.5/", "SIMULATED_0.5_figs", title="SIMULATED SIMPLE dataset, 0.5 errors per read", show_coords=True)
    plotplot(igrec_dir + "py/test_on_pd/SYNTHETIC_1/", "SYNTHETIC_1_figs", title="SYNTHETIC SIMPLE dataset, 1 error per read", show_coords=True)
    plotplot(igrec_dir + "py/test_on_pd/SYNTHETIC_2/", "SYNTHETIC_2_figs", title="SYNTHETIC SIMPLE dataset, 2 errors per read", show_coords=True)
    plotplot(igrec_dir + "py/test_on_pd/SYNTHETIC_0.5/", "SYNTHETIC_0.5_figs", title="SYNTHETIC SIMPLE dataset, 0.5 errors per read", show_coords=True)

    plotplot(igrec_dir + "py/test_on_pd/SIMULATED_1/", "Fig_11_a", title="SIMULATED SIMPLE dataset, 1.0 errors per read", show_coords=True)
    plotplot(igrec_dir + "py/test_on_pd/SIMULATED_0.5/", "Fig_11_b", title="SIMULATED SIMPLE dataset, 0.5 errors per read", show_coords=True)
    plotplot(igrec_dir + "py/test_on_pd/REAL/", "Fig_13", title="Sensitivity-precision plot (REAL dataset)", show_coords=True)
    # plotplot(igrec_dir + "py/flu_all/", "FLU_FV_21_IGH_3", title="")
