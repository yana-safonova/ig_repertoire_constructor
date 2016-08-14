from Bio import SeqIO
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

CUTAIL = 0

# TODO fix fonts
# FROM http://stackoverflow.com/questions/11367736/matplotlib-consistent-font-using-latex


current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = current_dir + "/../"
sys.path.append(igrec_dir + "/src/ig_tools/python_utils")
sys.path.append(igrec_dir + "/src/python_pipeline/")
import support
sys.path.append(igrec_dir + "/src/extra/ash_python_utils/")
from ash_python_utils import idFormatByFileName, smart_open, mkdir_p, FakeLog, memoize, memoize_invalidate

sys.path.append(igrec_dir + "/py")
from ig_compress_equal_clusters import parse_cluster_mult


path_to_igrec = igrec_dir


def hexbinblue(ax, x, y, gridsize=25, xlim=None, ylim=None):
    import matplotlib.colors as colors
    import matplotlib.pyplot as plt
    h = ax.hexbin(x, y,
                  gridsize=gridsize,
                  # bins='log',
                  norm=colors.LogNorm(),
                  mincnt=1,
                  edgecolor="darkblue",
                  linewidths=1,
                  cmap="Blues")
    plt.colorbar(h, ax=ax)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    xsize = (xlim[1] - xlim[0]) / float(gridsize)
    ysize = (ylim[1] - ylim[0]) / float(gridsize)
    ax.set_xlim((xlim[0] - xsize / 1.5), xlim[1] + xsize / 1.5)
    ax.set_ylim((ylim[0] - ysize / 1.5), ylim[1] + ysize / 1.5)

    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.set_aspect(abs(x1 - x0) / abs(y1 - y0))


def run_ig_matcher2(reference_file, constructed_file, output_file, prefix="", log=None,
                    tau=4, k=10, strategy=3):
    if log is None:
        log = FakeLog()

    if tau >= 9999:
        tau = 9999
        strategy = 0

    args = {"path": path_to_igrec,
            "swg_cmd": path_to_igrec + "/build/release/bin/ig_swgraph_construct",
            "reference_file": reference_file,
            "constructed_file": constructed_file,
            "output_file": output_file,
            "prefix": prefix,
            "k": k,
            "tau": tau,
            "strategy": strategy}

    support.sys_call("%(swg_cmd)s -r %(reference_file)s -i %(constructed_file)s -o %(output_file)s -k %(k)d --tau %(tau)d -A --strategy=%(strategy)s" % args,
                     log=log)


def run_consensus_constructor(initial_reads, rcm_file, output_file, log=None):
    if log is None:
        log = FakeLog()

    args = {"path": path_to_igrec,
            "cc_cmd": path_to_igrec + '/build/release/bin/ig_consensus_finder',
            "initial_reads": initial_reads,
            "output_file": output_file,
            "rcm_file": rcm_file}

    support.sys_call("%(cc_cmd)s -i %(initial_reads)s -R %(rcm_file)s -o %(output_file)s -H " % args,
                     log=log)


class MultToMultData:

    def __init__(self, constructed_abundances, constructed_sum, tau=0, reversed=False):
        import numpy as np

        constructed2reference = zip(constructed_abundances, constructed_sum[:, tau]) if not reversed else zip(constructed_sum[:, tau], constructed_abundances)
        constructed2reference = [x for x in constructed2reference if min(x[0], x[1]) > 0]
        constructed2reference.sort(key=lambda x: x[1])
        # TODO Proper fix
        __tmp = map(lambda x: np.array(x, dtype=float), zip(*constructed2reference))
        if len(__tmp) > 0:
            constructed_cluster_sizes, reference_cluster_sizes = __tmp
        else:
            constructed_cluster_sizes, reference_cluster_sizes = np.array([]), np.array([])

        rates = constructed_cluster_sizes / reference_cluster_sizes

        from cumulative_median import reversed_cumulative_median, reversed_cumulative_mean, unique
        median_rates = reversed_cumulative_median(rates)
        mean_rates = reversed_cumulative_mean(rates)

        # make theoretical abundances unique
        mean_rates_unique = unique(reference_cluster_sizes, mean_rates)
        median_rates_unique = unique(reference_cluster_sizes, median_rates)
        reference_cluster_sizes_unique = unique(reference_cluster_sizes)

        self.mean_rates_unique = np.array(mean_rates_unique)
        self.median_rates_unique = np.array(median_rates_unique)
        self.reference_cluster_sizes_unique = np.array(reference_cluster_sizes_unique)

        self.reference_cluster_sizes = reference_cluster_sizes
        self.constructed_cluster_sizes = constructed_cluster_sizes

    def median_rate(self, size=1):
        from bisect import bisect_left
        if len(self.reference_cluster_sizes_unique) == 0:
            return 0.

        i = bisect_left(self.reference_cluster_sizes_unique, size)  # the leftmost elt >= size
        assert i < len(self.reference_cluster_sizes_unique)
        return self.median_rates_unique[i]

    def mean_rate(self, size=1):
        from bisect import bisect_left
        if len(self.reference_cluster_sizes_unique) == 0:
            return 0.

        i = bisect_left(self.reference_cluster_sizes_unique, size)  # the leftmost elt >= size
        assert i < len(self.reference_cluster_sizes_unique)
        return self.mean_rates_unique[i]

    def plot_reference_vs_constructed_size(self, out, title="", format=None,
                                           points=True,
                                           marginals=False):
        import seaborn as sns

        f, ax = initialize_plot()

        if len(self.reference_cluster_sizes) == 0:
            # TODO Process more gently
            print "Could not plot"
            return None

        uplimit = max(self.reference_cluster_sizes)

        def round_up(number, ndigits=0):
            from math import ceil
            p = 10 ** ndigits
            return ceil(number * p) / p
        uplimit = round_up(uplimit, -2)

        if marginals:
            g = sns.JointGrid(x=self.reference_cluster_sizes,
                              y=self.constructed_cluster_sizes,
                              xlim=(0, uplimit), ylim=(0, uplimit),
                              ratio=5).set_axis_labels("Reference cluster size", "Constructed cluster size")

            g = g.plot_joint(sns.plt.scatter, alpha=0.5, label="clusters")
            g.plot_marginals(sns.distplot, kde=True, color=".5")

            ax = g.ax_joint
        else:
            if points:
                plt.plot(self.reference_cluster_sizes, self.constructed_cluster_sizes, "bo",
                         label="clusters", alpha=0.5,
                         markersize=6)

            else:
                f, ax = initialize_plot(figsize=(7.5, 6))
                # ax.set_axis_bgcolor('darkgrey')
                hexbinblue(ax, self.reference_cluster_sizes, self.constructed_cluster_sizes,
                           xlim=(0, uplimit), ylim=(0, uplimit))
                # plt.axes().set_aspect('equal', 'datalim')
                plt.axes().set_aspect('equal')

            ax.set_xlabel("Reference cluster size")
            ax.set_ylabel("Constructed cluster size")

        ax.plot([0, uplimit], [0, uplimit], "--", linewidth=0.5, color="black", label="$y = x$")
        # ax.plot([0, uplimit], [0, uplimit * self.median_rate(5)])

        # ax.plot(self.reference_cluster_sizes_unique, self.reference_cluster_sizes_unique * self.mean_rates_unique)
        ax.plot(self.reference_cluster_sizes_unique, self.reference_cluster_sizes_unique * self.median_rates_unique,
                label="rate median smoothing")

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc=2)

        if title:
            plt.title(title)

        save_plot(out, format=format)


class Reperoire2RepertoireMatching:

    def __init__(self,
                 constructed_repertoire, reference_repertoire,
                 tmp_file=None, max_tau=4, log=None):
        if tmp_file is None:
            import tempfile
            tmp_file = tempfile.mkstemp(suffix=".graph", prefix="aimquast_")[1]

        run_ig_matcher2(reference_repertoire, constructed_repertoire,
                        tau=max_tau,
                        output_file=tmp_file, log=log, prefix="")
        self.reference_abundances = get_clusters_sizes(reference_repertoire)

        ref_len = len(self.reference_abundances)

        self.constructed_abundances = []

        with smart_open(tmp_file) as f:
            header = f.next()
            cons_len, E, FORMAT = map(int, header.strip().split())

            self.constructed2reference = [[] for _ in xrange(cons_len)]
            self.reference2constructed = [[] for _ in xrange(ref_len)]

            for i, l in enumerate(f):
                l = l.strip().split()
                l = map(int, l)
                constructed_abundance = l[0]
                self.constructed_abundances.append(constructed_abundance)
                l = l[1:]
                neibs = [n - 1 for n in l[0::2]]
                dists = l[1::2]

                for j, d in zip(neibs, dists):
                    self.constructed2reference[i].append((j, d))
                    self.reference2constructed[j].append((i, d))

        import os
        os.remove(tmp_file)

    @staticmethod
    def bidirection(constructed_repertoire, reference_repertoire,
                    tmp_file=None, max_tau=4, log=None):
        r2c = Reperoire2RepertoireMatching(constructed_repertoire=constructed_repertoire,
                                           reference_repertoire=reference_repertoire,
                                           tmp_file=tmp_file, max_tau=max_tau, log=log)
        c2r = Reperoire2RepertoireMatching(constructed_repertoire=reference_repertoire,
                                           reference_repertoire=constructed_repertoire,
                                           tmp_file=tmp_file, max_tau=max_tau, log=log)

        def merge(x, y):
            assert len(x) == len(y)
            for i in xrange(len(x)):
                x[i] = list(set(x[i] + y[i]))  # TODO merge distances properly, use min

        merge(r2c.constructed2reference, c2r.reference2constructed)
        merge(r2c.reference2constructed, c2r.constructed2reference)

        return r2c

    def check(self,
              constructed_repertoire,
              reference_repertoire,
              log=None):
        if log is None:
            log = FakeLog()

        cons = []
        for i, neibs in enumerate(self.constructed2reference):
            matches = [j for j, d in neibs if d == 0]
            if len(matches) > 1:
                log.info("Cons %d matched on several references: %s" % (i, str(matches)))
                cons.append((i, matches))

        ref = []
        for j, neibs in enumerate(self.reference2constructed):
            matches = [i for i, d in neibs if d == 0]
            if len(matches) > 1:
                log.info("Ref %d matched on several constructed: %s" % (j, str(matches)))
                # print reference_repertoire[j]
                for i in matches:
                    # print constructed_repertoire[i]
                ref.append((j, matches))

        return cons, ref

    def plot_multiplicity_distributions(self,
                                        out,
                                        bins=25,
                                        ylog=True,
                                        format=None):
        import numpy as np
        import matplotlib.pyplot as plt

        f, ax = initialize_plot()

        _, bins = np.histogram(self.constructed_abundances + self.reference_abundances,
                               bins=bins)

        constructed_h, _ = np.histogram(self.constructed_abundances,
                                        bins=bins)
        reference_h, _ = np.histogram(self.reference_abundances,
                                      bins=bins)

        width = (bins[1] - bins[0]) / 3
        delta_outer = -width / 7
        delta_inner = width / 10

        if ylog:
            plt.yscale("log", nonposy="clip")

        ax.bar(bins[:-1] - width + delta_outer,
               reference_h,
               width=width - delta_outer - delta_inner,
               facecolor='cornflowerblue',
               label="Reference abundancies")

        ax.bar(bins[:-1] + delta_inner,
               constructed_h,
               width=width - delta_outer - delta_inner,
               facecolor='seagreen',
               label="Constructed abundancies")

        # plt.xticks(range(max_val + 1), labels)
        xlim = ax.get_xlim()
        xlim = (-width, xlim[1])
        ax.set_xlim(xlim)

        handles, labels = ax.get_legend_handles_labels()
        plt.legend(handles, labels)

        save_plot(out, format=format)


class RepertoireMatch:

    @staticmethod
    def __read_fa(filename):
        with smart_open(filename) as f:
            result = map(lambda rec: str(rec.seq), SeqIO.parse(f, idFormatByFileName(filename)))

        return result

    @staticmethod
    def __errs(s1, s2):
        result = []
        for i, l1, l2 in zip(xrange(1000), s1, s2):
            if l1 != l2:
                result.append(i)

        return result

    @staticmethod
    def __hamming(s1, s2):
        result = 0
        for l1, l2 in zip(s1, s2):
            if l1 != l2:
                result += 1

        return result

    @memoize
    def __error_positions(self, tau=1, size=5,
                          reversed=True, relative=False):
        result = []

        for i in xrange(len(self.reference_repertoire)):
            if self.reference[i, tau - 1] < size and self.reference[i, tau] >= size:
                j = self.ref_edge[i, tau]
                l = len(self.reference_repertoire[i])

                if j > -1:
                    positions = self.__errs(self.reference_repertoire[i], self.constructed_repertoire[j])
                    if reversed:
                        positions = map(lambda pos: l - pos, positions)
                    if relative:
                        positions = map(lambda pos: float(pos) / float(l), positions)
                    result.extend(positions)

        result.sort()
        return result


    def __init__(self,
                 constructed_repertoire, reference_repertoire,
                 tmp_file=None, max_tau=10,
                 reference_trash_cutoff=-float("inf"),
                 reference_trust_cutoff=float("inf"),
                 log=None):
        import numpy as np

        self.reference_repertoire = self.__read_fa(reference_repertoire)
        self.constructed_repertoire = self.__read_fa(constructed_repertoire)

        self.rep2rep = Reperoire2RepertoireMatching.bidirection(constructed_repertoire=constructed_repertoire,
                                                                reference_repertoire=reference_repertoire,
                                                                max_tau=max_tau,
                                                                tmp_file=tmp_file,
                                                                log=log)

        self.rep2rep.check(log=log,
                           reference_repertoire=self.reference_repertoire,
                           constructed_repertoire=self.constructed_repertoire)

        assert reference_trash_cutoff <= reference_trust_cutoff

        ref_len = len(self.rep2rep.reference_abundances)
        cons_len = len(self.rep2rep.constructed_abundances)

        reference = np.full((ref_len, max_tau + 1), 0, dtype=int)
        constructed = np.full((cons_len, max_tau + 1), 0, dtype=int)

        ref_edge = np.full((ref_len, max_tau + 1), -1, dtype=int)
        cons_edge = np.full((cons_len, max_tau + 1), -1, dtype=int)

        reference_sum = np.full((ref_len, max_tau + 1), 0, dtype=int)
        constructed_sum = np.full((cons_len, max_tau + 1), 0, dtype=int)

        assert len(self.rep2rep.constructed2reference) == len(self.constructed_repertoire)
        assert len(self.rep2rep.reference2constructed) == len(self.reference_repertoire)

        for i, pairs in enumerate(self.rep2rep.constructed2reference):
            for j, d in pairs:
                assert self.__hamming(self.constructed_repertoire[i], self.reference_repertoire[j]) == d

                reference_abundance = self.rep2rep.reference_abundances[j]
                constructed_abundance = self.rep2rep.constructed_abundances[i]

                reference_sum[j, d] = reference_sum[j, d] + constructed_abundance
                constructed_sum[i, d] = constructed_sum[i, d] + reference_abundance

                if reference_abundance >= reference_trust_cutoff:
                    reference_abundance = float("inf")
                if reference_abundance < reference_trash_cutoff:
                    reference_abundance = -float("inf")

                min_abundance = min(constructed_abundance, reference_abundance)
                if reference[j, d] < min_abundance:
                    reference[j, d] = min_abundance
                    ref_edge[j, d] = i

                if constructed[i, d] < min_abundance:
                    constructed[i, d] = min_abundance
                    cons_edge[i, d] = j

        for d in xrange(1, max_tau + 1):
            for j in xrange(ref_len):
                reference[j, d] = max(reference[j, d], reference[j, d - 1])
                reference_sum[j, d] = reference_sum[j, d] + reference_sum[j, d - 1]

            for i in xrange(cons_len):
                constructed[i, d] = max(constructed[i, d], constructed[i, d - 1])
                constructed_sum[i, d] = constructed_sum[i, d] + constructed_sum[i, d - 1]

        self.cons_edge = cons_edge
        self.ref_edge = ref_edge
        self.constructed = constructed
        self.reference = reference

        self.M2MDATA = MultToMultData(self.rep2rep.reference_abundances, reference_sum, reversed=True)

        self.sensitivity_vectors = []
        for col in reference.T:
            self.sensitivity_vectors.append(sorted(col))

        self.precision_vectors = []
        for col in constructed.T:
            self.precision_vectors.append(sorted(col))

        self.constructed_abundances = sorted(self.rep2rep.constructed_abundances)
        self.reference_abundances = sorted(self.rep2rep.reference_abundances)

        self.max_tau = max_tau
        self.reference_trust_cutoff = reference_trust_cutoff
        self.reference_trash_cutoff = reference_trash_cutoff

    def plot_reference_vs_constructed_size(self, *args, **kwargs):
        return self.M2MDATA.plot_reference_vs_constructed_size(*args, **kwargs)

    def plot_multiplicity_distributions(self, *args, **kwargs):
        return self.rep2rep.plot_multiplicity_distributions(*args, **kwargs)

    def reference_size(self, size=1):
        from cumulative_median import how_many_greater_or_equal

        size = min(size, self.reference_trust_cutoff)
        size = max(size, self.reference_trash_cutoff)

        assert size > 0
        return how_many_greater_or_equal(self.reference_abundances, size)

    def constructed_size(self, size=1):
        from cumulative_median import how_many_greater_or_equal

        assert size > 0
        return how_many_greater_or_equal(self.constructed_abundances, size)

    def ref2cons(self, size=1, tau=0):
        from cumulative_median import how_many_greater_or_equal

        assert size > 0
        assert 0 <= tau <= self.max_tau

        return how_many_greater_or_equal(self.sensitivity_vectors[tau], size)

    def sensitivity(self, size=1, tau=0):
        all = self.reference_size(size)
        identified = self.ref2cons(size, tau)

        assert(all >= identified)
        return float(identified) / float(all) if all > 0 else 0.

    def cons2ref(self, size=1, tau=0):
        from cumulative_median import how_many_greater_or_equal

        assert size > 0
        assert 0 <= tau <= self.max_tau

        return how_many_greater_or_equal(self.precision_vectors[tau], size)

    def precision(self, size=1, tau=0):
        all = self.constructed_size(size)
        true = self.cons2ref(size, tau)

        assert(all >= true)
        return float(true) / float(all) if all > 0 else 0.

    def fdr(self, size=1, tau=0):
        return 1 - self.precision(size, tau)

    def f1(self, size=1, tau=0):
        p = self.precision(size, tau)
        r = self.sensitivity(size, tau)
        return 2 * (p * r) / (p + r)

    # TODO Add density plot(s)

    def report(self, report, size=None):
        if "reference_based" not in report.__dict__:
            report.reference_based = {}

        rb = report.reference_based

        if size is None:
            size = self.reference_trust_cutoff

        rb["min_size"] = size

        rb["precision"] = self.precision(size)
        rb["sensitivity"] = self.sensitivity(size)
        rb["cons2ref"] = self.cons2ref(size)
        rb["ref2cons"] = self.ref2cons(size)

        rb["constructed_size"] = self.constructed_size(size)
        rb["reference_size"] = self.reference_size(size)

        rb["reference_vs_constructed_size_median_rate"] = float(self.M2MDATA.median_rate(size))
        rb["reference_vs_constructed_size_mean_rate"] = float(self.M2MDATA.mean_rate(size))  # TODO check type

        precision, sizes = self.__get_data(what="precision")
        sensitivity, _ = self.__get_data(what="sensitivity")
        rb["__data_sizes"] = list(sizes)
        rb["__data_sensitivity"] = list(sensitivity)
        rb["__data_precision"] = list(precision)

        ref2cons5, taus, _ = self.__get_measure_for_plotting(size=5, what="ref2cons")
        cons2ref5, _, _ = self.__get_measure_for_plotting(size=5, what="cons2ref")
        rb["__data_taus"] = list(taus)
        rb["__data_ref2cons5"] = list(ref2cons5)
        rb["__data_cons2ref5"] = list(cons2ref5)

        rb["error_positions_reversed"] = self.__error_positions(reversed=True)
        rb["error_positions"] = self.__error_positions(reversed=False)

    def __get_measure_for_plotting(self,
                                   size=1,
                                   what="sensitivity",
                                   differential=True,
                                   max_tau=4):
        assert what in ["sensitivity", "fdr", "precision", "recall", "ref2cons", "cons2ref"]

        max_tau = min(max_tau, self.max_tau)

        def func(tau):
            return getattr(self, what)(size=size, tau=tau)

        taus = range(max_tau + 1)
        measures = map(func, taus)
        labels = ["$%d$" % i for i in taus]

        if what == "ref2cons":
            ALL = self.reference_size(size)
        elif what == "cons2ref":
            ALL = self.constructed_size(size)
        else:
            ALL = 1

        if differential:
            if what != "fdr":
                residual = ALL - measures[len(measures) - 1]
                for i in reversed(xrange(1, len(measures))):
                    measures[i] -= measures[i - 1]
                measures.append(residual)
                labels.append(r"$\geq %d$" % (max_tau + 1))
                taus.append(max_tau + 1)
            else:
                for i in xrange(1, len(measures)):
                    measures[i - 1] -= measures[i]
                labels[len(measures) - 1] = r"$\geq %d$" % (max_tau)

        return measures, taus, labels

    def plot_sensitivity_precision(self, out,
                                   size=1,
                                   what="sensitivity",
                                   differential=True,
                                   max_tau=None,
                                   format=None):
        import numpy as np
        import matplotlib.pyplot as plt

        f, ax = initialize_plot(font_scale=1)

        if max_tau is None:
            max_tau = self.max_tau

        data = self.__get_measure_for_plotting(size=size, what=what, differential=differential, max_tau=max_tau)
        measures = np.array(data[0])
        taus = np.array(data[1])
        labels = data[2]

        width = 0.9
        ax.bar(taus + 0.5 - width / 2, measures, width=width,
               facecolor='cornflowerblue',
               label="Actual frequencies")
        plt.xticks(taus + 0.5,
                   labels)

        plt.xlim(0, max(taus) + 1)
        plt.ylim(0, sum(measures))

        if what in ["sensitivity", "ref2cons"]:
            plt.title("Distribution of distance from reference to constructed sequences")
        elif what in ["precision", "cons2ref"]:
            plt.title("Distribution of distance from constructed to reference sequences")

        ax.set_xlabel("distance")
        ax.set_ylabel("#clusters")
        if what not in ["ref2cons", "cons2ref"]:
            plt.ylim(0, 1)

        save_plot(out, format=format)

    def plot_octoplot(self, out,
                      sizes=(1, 3, 5, 10),
                      differential=True,
                      max_tau=None,
                      format=None):
        import numpy as np
        import matplotlib.pyplot as plt

        if max_tau is None:
            max_tau = self.max_tau

        initialize_plot(figsize=(8, 4), font_scale=1)

        N = len(sizes)
        for i in xrange(len(sizes)):
            size = sizes[i]
            plt.subplot(2, N, i + 1)

            data = self.__get_measure_for_plotting(size=size,
                                                   what="ref2cons",
                                                   differential=differential,
                                                   max_tau=max_tau)
            measures = np.array(data[0])
            taus = np.array(data[1])
            labels = data[2]

            width = 0.9
            plt.bar(taus + 0.5 - width / 2, measures, width=width,
                    facecolor='cornflowerblue')
            plt.xticks(taus + 0.5,
                       labels)
            plt.title("sensitivity, size = %d" % size)
            plt.xlim(0, max(taus) + 1)
            plt.ylim(0, sum(measures))

            plt.subplot(2, N, N + i + 1)

            data = self.__get_measure_for_plotting(size=size,
                                                   what="cons2ref",
                                                   differential=differential,
                                                   max_tau=max_tau)
            measures = np.array(data[0])
            taus = np.array(data[1])
            labels = data[2]

            width = 0.9
            plt.bar(taus + 0.5 - width / 2, measures, width=width,
                    facecolor='orangered')
            plt.xticks(taus + 0.5,
                       labels)
            plt.title("precision, size = %d" % size)
            plt.xlim(0, max(taus) + 1)
            plt.ylim(0, sum(measures))

        plt.tight_layout()

        save_plot(out, format=format)

    def __get_data(self,
                   what="precision",
                   max_size_threshold=75,
                   tau=0):
        import numpy as np

        assert tau <= self.max_tau

        assert what in ["sensitivity", "fdr", "precision", "recall"]

        def func(size):
            return getattr(self, what)(size=size, tau=tau)

        sizes = np.array(range(1, max_size_threshold + 1))

        x = np.array(map(func, sizes))

        return x, sizes

    def plot_error_pos_dist(self, out, format=None):
        import seaborn as sns

        f, ax = initialize_plot()

        errors = self.__error_positions(reversed=False, relative=True)
        try:
            sns.distplot(errors, kde=False, bins=25, ax=ax)
            ax.set_xlabel("Relative error positions")
            ax.set_xlim((0., 1.))

            save_plot(out, format=format)
        except BaseException as ex:
            print ex

    def plot_min_cluster_size_choose(self,
                                     what_x="precision",
                                     what_y="sensitivity",
                                     max_size_threshold=75,
                                     tau=0,
                                     out="fdr_sensitivity",
                                     title="",
                                     format=None):
        import numpy as np
        import matplotlib.pyplot as plt
        import seaborn as sns

        f, ax = initialize_plot()

        x, sizes = self.__get_data(what=what_x, max_size_threshold=max_size_threshold, tau=tau)
        y, _ = self.__get_data(what=what_y, max_size_threshold=max_size_threshold, tau=tau)

        sns.set_style("darkgrid")
        plt.plot(x, y, "b-", color="cornflowerblue")

        eps = 0.025
        plt.xlim((0 - eps, 1 + eps))
        plt.ylim((0 - eps, 1 + eps))

        plt.ylabel(what_y)
        plt.xlabel(what_x)

        blue_points = []

        def annotation(size, text=None, xshift=-0, yshift=-0.05):
            if text is None:
                text = str(size)

            _x = x[size - 1]
            _y = y[size - 1]

            bp = plt.plot(_x, _y, "bo", color="blue", label="min size")
            blue_points.append(bp)

            plt.annotate(text, xy=(_x, _y),
                         xytext=(_x + xshift, _y + yshift))

        annotation(1)
        annotation(3)
        annotation(5)
        annotation(10)
        annotation(50)

        reference_trust_cutoff = self.reference_trust_cutoff
        if np.isfinite(reference_trust_cutoff):
            reference_trust_cutoff = int(reference_trust_cutoff)
            mask = sizes == reference_trust_cutoff
            plt.plot(x[mask], y[mask], "bo", color="red", label="Reference size threshold")

        if title:
            plt.title(title)

        handles, labels = ax.get_legend_handles_labels()
        handles = [handles[0], handles[-1]]
        labels = [labels[0], labels[-1]]
        ax.legend(handles, labels, loc=3)

        save_plot(out, format=format)


def get_clusters_sizes(filename):
    sizes = []
    with smart_open(filename) as fin:
        for record in SeqIO.parse(fin, idFormatByFileName(filename)):
            clm = parse_cluster_mult(record.description)
            if clm is None:
                cluster, mult = None, 1
            else:
                cluster, mult = parse_cluster_mult(record.description)
            sizes.append(mult)
    return sizes


def parse_rcm(filename):
    rcm = {}
    with smart_open(filename) as f:
        for line in f:
            line = [_.strip() for _ in line.strip().split("\t")]
            if len(line) == 2:
                id, cluster = line
            else:
                id, cluster = line[0], None

            rcm[id] = cluster

    # TODO compute cluster size
    return rcm


class RcmVsRcm:

    @staticmethod
    def rcm_count_cluster_sizes(rcm):
        from collections import defaultdict

        cluster_sizes = defaultdict(int)
        for cluster in rcm.itervalues():
            cluster_sizes[cluster] += 1

        return cluster_sizes

    def __init__(self, rcm1, rcm2):
        rcm1, rcm2 = parse_rcm(rcm1), parse_rcm(rcm2)

        self.sizes1 = self.rcm_count_cluster_sizes(rcm1)
        self.sizes2 = self.rcm_count_cluster_sizes(rcm2)

        ids = list(set(rcm1.keys() + rcm2.keys()))

        self.clustering1 = []
        self.clustering2 = []

        for id in ids:
            self.clustering1.append(rcm1[id] if id in rcm1 else None)
            self.clustering2.append(rcm2[id] if id in rcm2 else None)

    @memoize
    def indices(self):
        from copy import copy

        def fix_nones(v, prefix):
            for i in xrange(len(v)):
                if v[i] is None:
                    v[i] = prefix + str(i)

        clustering1 = copy(self.clustering1)
        clustering2 = copy(self.clustering2)
        fix_nones(clustering1, "__none__clustering1__")
        fix_nones(clustering2, "__none__clustering2__")

        indices = clustering_similarity_indices(clustering1, clustering2)
        indices.constructed_purity = purity(clustering1, clustering2)
        indices.reference_purity = purity(clustering2, clustering1)

        return indices

    @memoize_invalidate
    def prune(self, size1=1, size2=1):
        N = len(self.clustering1)

        for i in xrange(N):
            if self.sizes1[self.clustering1[i]] < size1:
                self.clustering1[i] = None
            if self.sizes2[self.clustering2[i]] < size2:
                self.clustering2[i] = None

    def prune_copy(self, size1=1, size2=1):
        from copy import deepcopy
        result = deepcopy(self)
        result.prune(size1, size2)
        return result

    def report(self, report, name="reference_based"):
        if name not in report.__dict__:
            report.__dict__[name] = {}

        rb = report.__dict__[name]

        indices = self.indices()
        for k, v in indices.__dict__.iteritems():
            rb[k] = v

        constructed_purity = self.purity(constructed=True)
        reference_purity = self.purity(constructed=False)
        rb["n_dirty_constructed_clusters"] = sum(constructed_purity < 0.95)
        rb["n_dirty_reference_clusters"] = sum(reference_purity < 0.95)
        rb["n_pure_constructed_clusters"] = sum(constructed_purity >= 0.95)
        rb["n_pure_reference_clusters"] = sum(reference_purity >= 0.95)

        rb["constructed_secondary_votes"] = sum(self.secondary_votes(constructed=True))
        rb["reference_secondary_votes"] = sum(self.secondary_votes(constructed=False))

    @memoize
    def votes(self, constructed=True):
        return votes(self.clustering1, self.clustering2) if constructed else votes(self.clustering2, self.clustering1)

    @memoize
    def purity(self, constructed=True):
        import numpy as np

        votes = self.votes(constructed)
        majority_votes = [vote[0] for vote in votes]
        # secondary_votes = [vote[1] for vote in votes]
        sizes = [sum(vote) for vote in votes]

        return np.array(majority_votes, dtype=float) / np.array(sizes, dtype=float)

    @memoize
    def majority_votes(self, constructed=True):
        votes = self.votes(constructed)
        majority_votes = [vote[0] for vote in votes]

        return majority_votes

    @memoize
    def secondary_votes(self, constructed=True):
        votes = self.votes(constructed)
        secondary_votes = [vote[1] for vote in votes]

        return secondary_votes

    def plot_majority_secondary(self, out, format=None, constructed=True,
                                points=False):
        import numpy as np
        f, ax = initialize_plot()

        votes = self.votes(constructed)
        majority_votes = np.array([vote[0] for vote in votes])
        secondary_votes = np.array([vote[1] for vote in votes])

        try:
            if points:
                ax.plot(majority_votes, secondary_votes, "bo",
                        alpha=0.4,
                        markersize=6,
                        label="clusters")
            else:
                # ax.set_axis_bgcolor('darkgrey')
                f, ax = initialize_plot(figsize=(7.5, 6))
                hexbinblue(ax, majority_votes, secondary_votes)
                # cb = plt.colorbar()
                # cb.set_label('log10(N)')
            ax.set_xlabel("Majority votes")  # Primary
            ax.set_ylabel("Secondary votes")
        except BaseException as ex:
            print ex

        save_plot(out, format=format)

    def plot_size_nomajority(self, out, format=None, constructed=True,
                             points=False):
        import numpy as np

        f, ax = initialize_plot()

        votes = self.votes(constructed)
        majority_votes = np.array([vote[0] for vote in votes])
        sizes = np.array([sum(vote) for vote in votes])

        try:
            # sns.kdeplot(sizes, sizes - majority_votes, ax=ax,
            #             cmap="Blues",
            #             shade_lowest=False,
            #             shade=True)
            if points:
                ax.plot(sizes, sizes - majority_votes, "bo",
                        alpha=0.4,
                        markersize=1,
                        label="clusters")
            else:
                # y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
                # ax.yaxis.set_major_formatter(y_formatter)
                # ax.xaxis.set_major_formatter(y_formatter)
                f, ax = initialize_plot(figsize=(7.5, 6))
                hexbinblue(ax, sizes, sizes - majority_votes)
            ax.set_xlabel("Cluster size")  # Primary
            ax.set_ylabel("Cluster size - majority votes")

            save_plot(out, format=format)
        except BaseException as ex:
            print ex

    def plot_purity_distribution(self, out, format=None, constructed=True):
        import seaborn as sns

        f, ax = initialize_plot()

        purity = self.purity(constructed)
        try:
            sns.distplot(purity, kde=False, bins=25, ax=ax)
            ax.set_ylabel("Purity")
            ax.set_xlim((0, 1))
            ax.set_ylim((0, len(purity)))

            save_plot(out, format=format)
        except BaseException as ex:
            print ex


def reconstruct_rcm(initial_reads, repertoire,
                    tmp_file_matcher=None, tmp_file_reads=None,
                    taus=(1, 2, 4, 8, 12),
                    fallback_to_exhaustive_mode=False,
                    log=None):
    if log is None:
        log = FakeLog()

    import random

    from collections import defaultdict
    rcm = defaultdict(str)

    # get cluster names
    with smart_open(repertoire) as f:
        cluster_names = [parse_cluster_mult(record.id)[0] for record in SeqIO.parse(f, idFormatByFileName(repertoire))]

    # get read ids
    with smart_open(initial_reads) as f:
        read_ids = [str(record.description) for record in SeqIO.parse(f, idFormatByFileName(initial_reads))]

    unmatched_reads = set(read_ids)
    if tmp_file_reads is None:
        import tempfile
        tmp_file_reads = tempfile.mkstemp(suffix=".fa", prefix="input_reads_")[1]

    uncertain = 0

    taus = list(set(taus))
    taus.sort()
    if fallback_to_exhaustive_mode:
        taus.append(float("inf"))

    for tau in taus:
        # Save unmatched reads to tmp file
        written_reads = 0
        read_ids = []
        with smart_open(initial_reads) as fin, smart_open(tmp_file_reads, "w") as fout:
            for record in SeqIO.parse(fin, idFormatByFileName(initial_reads)):
                if str(record.description) in unmatched_reads:
                    SeqIO.write(record, fout, idFormatByFileName(tmp_file_reads))
                    written_reads += 1
                    read_ids.append(record.description)

        log.info("%d reads are written to tmp file" % written_reads)
        r2r = Reperoire2RepertoireMatching.bidirection(reference_repertoire=tmp_file_reads,
                                                       tmp_file=tmp_file_matcher,
                                                       constructed_repertoire=repertoire,
                                                       log=log, max_tau=tau)

        for read_index, l in enumerate(r2r.reference2constructed):
            if l:
                dists = [d for i, d in l]
                min_dist = min(dists)
                nearest_neibs = [i for i, d in l if d == min_dist]

                if len(nearest_neibs) > 1:
                    uncertain += 1

                assigned_neib = random.choice(nearest_neibs)
                rcm[read_ids[read_index]] = cluster_names[assigned_neib]
                unmatched_reads.remove(read_ids[read_index])

        log.info("%d unmatched reads left" % len(unmatched_reads))

        if not unmatched_reads:
            break

    log.info("Uncertain reads %d" % uncertain)

    import os
    os.remove(tmp_file_reads)
    return rcm


def write_rcm(rcm, filename):
    with smart_open(filename, "w") as f:
        for read, cluster in rcm.iteritems():
            f.write("%s\t%s\n" % (read, cluster))


def clustering_similarity_indices(X, Y):
    from collections import defaultdict
    from math import sqrt, log

    m = defaultdict(int)
    mX = defaultdict(int)
    mY = defaultdict(int)

    for x, y in zip(X, Y):
        m[(x, y)] += 1
        mX[x] += 1
        mY[y] += 1

    S00 = S01 = S10 = S11 = 0

    assert len(Y) == len(X)
    N = len(X)

    for x, y in zip(X, Y):
        _m = m[(x, y)]
        _mX = mX[x]
        _mY = mY[y]
        S00 += _m - 1
        S01 += _mX - _m
        S10 += _mY - _m
        S11 += N - _mX - _mY + _m

    S00 //= 2
    S01 //= 2
    S10 //= 2
    S11 //= 2

    def comb2(n):
        return n * (n - 1) // 2

    # Compute the ARI using the contingency data
    sum_comb_c = sum(comb2(_) for _ in mY.itervalues())
    sum_comb_k = sum(comb2(_) for _ in mX.itervalues())
    sum_comb = sum(comb2(_) for _ in m.itervalues())

    assert sum_comb == S00
    assert S01 == sum_comb_k - S00
    assert S10 == sum_comb_c - S00
    assert S11 == comb2(N) - S01 - S10 - S00

    prod_comb = (sum_comb_c * sum_comb_k) / float(comb2(N))
    mean_comb = (sum_comb_k + sum_comb_c) / 2.

    class SimIndices:

        def __init__(self):
            self.adjusted_rand_index = 1.
            self.rand_index = 1.
            self.fowlkes_mallows_index = 1.
            self.jaccard_index = 1.

    result = SimIndices()
    result.adjusted_rand_index = (sum_comb - prod_comb) / (mean_comb - prod_comb) if mean_comb - prod_comb > 0 else 1.
    result.rand_index = float(S00 + S11) / float(comb2(N))
    result.fowlkes_mallows_index = float(S00) / sqrt(float(S00 + S10) * float(S00 + S01)) if S00 + S10 > 0 and S00 + S01 else 1.
    result.jaccard_index = float(S00) / float(S00 + S10 + S01) if S00 + S10 + S01 else 1.

    # FROM http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
    mutual_information = 0.
    for key, count in m.iteritems():
        x, y = key
        mutual_information += float(count) / N * (log(N) + log(count) - log(mX[x]) - log(mY[y]))

    def entropy(freqs, N=None):
        if N is None:
            N = sum(freqs)
        return -sum(float(n) / float(N) * (log(n) - log(N)) for n in freqs)

    result.normalized_mutual_information = mutual_information / (entropy(mX.values()) + entropy(mY.values())) * 2.

    return result


def purity(X, Y):
    votes__ = votes(X, Y)
    majority_votes = sum(vote[0] for vote in votes__)

    return float(majority_votes) / float(len(X))


def votes(X, Y):
    from collections import defaultdict

    assert len(X) == len(Y)

    class defaultdict_factory:

        def __init__(self, type):
            self.__type = type

        def __call__(self):
            return defaultdict(self.__type)

    cluster = defaultdict(defaultdict_factory(int))  # Unfortunately, it's impossible to make defaultdict(defaultdict)
    for x, y in zip(X, Y):
        if x is not None and y is not None:
            cluster[x][y] += 1

    votes = [sorted(cluster_content.itervalues(), reverse=True) for cluster_content in cluster.itervalues()]

    for vote in votes:
        if len(vote) < 2:
            vote += [0] * (2 - len(vote))

    return votes


def dict2list(d, default=int):
    if not d:
        return []

    for k in d.iterkeys():
        assert type(k) is int
        assert k >= 0

    max_key = max(d.iterkeys())
    return [d[k] if k in d else default() for k in xrange(max_key + 1)]

assert dict2list({1: 10, 5: 7}) == [0, 10, 0, 0, 0, 7]


def initialize_plot(figsize=(6, 6), font_scale=1.5):
    import matplotlib.pyplot as plt
    import seaborn as sns
    plt.ioff()
    # matplotlib.use('Agg')

    # matplotlib.rcParams['mathtext.fontset'] = 'stix'
    # matplotlib.rcParams['font.family'] = 'STIXGeneral'

    # matplotlib.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    # matplotlib.rc('text', usetex=True)

    # matplotlib.rcParams['mathtext.fontset'] = 'custom'
    # matplotlib.rcParams['mathtext.tt'] = matplotlib.rcParams['font.monospace']
    # matplotlib.rcParams['mathtext.sf'] = matplotlib.rcParams['font.sans-serif']
    # matplotlib.rcParams['mathtext.rm'] = matplotlib.rcParams['font.serif']
    # # matplotlib.rcParams['mathtext.cal'] = matplotlib.rcParams[]
    # matplotlib.rcParams['mathtext.it'] = matplotlib.rcParams['font.cursive']
    # matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

    sns.set_style("darkgrid")
    sns.set(font_scale=font_scale)

    f, ax = plt.subplots(figsize=figsize)

    return f, ax


def save_plot(plot_name,
              format=None,
              make_dir=True,
              close=True):
    import matplotlib.pyplot as plt
    import os.path

    if format is None:
        format = ('png', 'pdf', 'svg')

    if make_dir:
        dirname = os.path.dirname(plot_name)
        if dirname:
            mkdir_p(dirname)

    plt.tight_layout()

    if 'png' in format:
        plt.savefig(plot_name + '.png', format='png')

    if 'pdf' in format:
        from matplotlib.backends.backend_pdf import PdfPages
        pp = PdfPages(plot_name + '.pdf')
        pp.savefig()
        pp.close()

    if 'svg' in format:
        plt.savefig(plot_name + '.svg', format='svg')

    if close:
        plt.close()


def consensus(reads):
    import numpy as np

    # l = max(len(read) for read in reads)
    l = min(len(read) for read in reads)

    nuc2idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    idx2nuc = np.array(list("ACGT"))

    reads = [str(read) for read in reads]

    mx = np.zeros(shape=(l, 4))
    for read in reads:
        for i in xrange(min(l, len(read))):
            letter = read[i]
            idx = nuc2idx[letter]
            mx[i, idx] += 1

    indices = mx.argmax(axis=1)
    consensus = idx2nuc[indices]
    consensus = "".join(consensus)

    # TODO Return misconsensus

    return consensus

assert consensus(["GAAA", "AAAC", "AATA"]) == "AAAA"


def majority_secondary(reads):
    import numpy as np

    # l = max(len(read) for read in reads)
    l = min(len(read) for read in reads)

    nuc2idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    reads = [str(read) for read in reads]

    mx = np.zeros(shape=(l, 4))
    for read in reads:
        for i in xrange(min(l, len(read))):
            letter = read[i]
            idx = nuc2idx[letter]
            mx[i, idx] += 1

    mx.sort(axis=1)
    return mx[:, 3], mx[:, 2]

assert consensus(["GAAA", "AAAC", "AATA"]) == "AAAA"


class Cluster:

    def __init__(self):
        self.__reads = []
        self.mult = self.center = self.name = None

    @memoize_invalidate
    def append(self, read):
        self.__reads.append(read)
        return self

    @memoize
    def consensus(self):
        reads = [str(read.seq) for read in self.__reads]
        return consensus(reads)

    def centroid(self):
        if self.center is not None:
            return self.center
        else:
            return self.consensus()

    def export(self, filename):
        with smart_open(filename, "w") as fout:
            for read in self.__reads:
                SeqIO.write(read, fout, idFormatByFileName(filename))

    def __len__(self):
        return len(self.__reads)

    @memoize
    def length(self):
        return min(len(read) for read in self.__reads)

    @memoize
    def majority_secondary(self):
        reads = [str(read.seq) for read in self.__reads]
        ms = majority_secondary(reads)
        l = self.length()
        return ms[0][:l], ms[1][:l]

    @memoize
    def nerrors_by_position(self):
        return self.majority_secondary()[1]

    @memoize
    def nerrors_by_read(self):
        reads = [str(read.seq) for read in self.__reads]
        center = self.centroid()
        errors_by_read = []

        for read in reads:
            error_by_read = 0

            for i in xrange(min(len(read), len(center) - CUTAIL)):
                if read[i] != center[i]:
                    error_by_read += 1

            errors_by_read.append(error_by_read)

        return errors_by_read

    @memoize
    def errors01(self):
        result = []
        errors_by_position = self.nerrors_by_position()

        l = self.length()
        for i, nerrors in enumerate(errors_by_position):
            if nerrors:
                result += [float(i) / float(l)] * nerrors

        return result

    @memoize
    def max_error(self):
        return max(self.nerrors_by_position()) if self.length() else 0

    @memoize
    def max_ppf(self, q, size, error_rate):
        from scipy.stats import poisson
        import math

        l = self.length()

        prob_error = float(error_rate) / float(l) / 3
        lam = float(prob_error) * float(size)
        return poisson.ppf(math.pow(q, 1. / l), mu=lam)

    @memoize
    def max_cdf(self, x, size, error_rate):
        from scipy.stats import poisson

        l = self.length()

        prob_error = float(error_rate) / float(l) / 3
        lam = float(prob_error) * float(size)
        return poisson.cdf(x, mu=lam) ** l

    @memoize
    def pvalue_upper(self, error_rate):
        return 1 - self.max_cdf(self.max_error() - 1, len(self), error_rate)

    @memoize
    def pvalue_lower(self, error_rate):
        F = self.max_cdf(self.max_error(), len(self), error_rate)
        return F

    @memoize
    def pvalue_both(self, error_rate):
        return 2 * min(self.pvalue_upper(error_rate), self.pvalue_lower(error_rate))

    @memoize
    def nerrors(self):
        return sum(self.nerrors_by_position())


def join_list_of_lists(a):
    import itertools
    return list(itertools.chain.from_iterable(a))


class Repertoire:

    def __init__(self, rcm, library, repertoire, min_size=5):
        from Bio import SeqIO
        from collections import defaultdict

        self.__min_size = min_size

        with smart_open(library) as f:
            reads = [rec for rec in SeqIO.parse(f, idFormatByFileName(library))]

        id2read = {str(rec.description): rec for rec in reads}

        rcm = parse_rcm(rcm)

        with smart_open(repertoire) as f:
            reads = [rec for rec in SeqIO.parse(f, idFormatByFileName(repertoire))]

        clusters = defaultdict(Cluster)

        for id, cluster in rcm.iteritems():
            if id in id2read:
                clusters[cluster].append(id2read[id])
                clusters[cluster].name = cluster

        for read in reads:
            cluster, mult = parse_cluster_mult(str(read.description))
            if cluster in clusters:
                clusters[cluster].mult = mult
                clusters[cluster].center = read
                clusters[cluster].name = cluster

        self.__clusters = clusters

    @memoize
    def error_rates(self, min_size=None):
        if min_size is None:
            min_size = self.__min_size

        nerrors = self.__nerrors_by_read(min_size)

        # Different strategies
        # MLE
        # first-len
        # By 1-2

        # Group
        from collections import defaultdict
        errordist = defaultdict(int)
        for e in nerrors:
            errordist[e] += 1

        N = len(nerrors)
        # MLE
        #

        class Empty:

            def __str__(self):
                return str(self.__dict__)

        res = Empty()
        res.mle = float(sum(nerrors)) / N

        # first - len MLE by ideal
        import math
        if errordist[0] > 0:
            res.first_len = - math.log(float(errordist[0]) / N)
        else:
            res.fist_len = res.mle

        # first_second
        if errordist[0] > 0:
            res.first_second = float(errordist[1]) / float(errordist[0])
        else:
            res.first_second = res.mle

        # first_third
        if errordist[0] > 0:
            res.first_third = math.sqrt(2. * float(errordist[2]) / float(errordist[0]))
        else:
            res.first_third = res.mle

        return res

    def error_rate(self, min_size=None):
        return self.error_rates(min_size=min_size).first_third

    @memoize
    def __nerrors_by_read(self, min_size=None):
        if min_size is None:
            min_size = self.__min_size

        return join_list_of_lists(cluster.nerrors_by_read() for cluster in self.__clusters.itervalues() if len(cluster) >= min_size)

    @memoize
    def __errors01(self, min_size=None):
        if min_size is None:
            min_size = self.__min_size

        return join_list_of_lists(cluster.errors01() for cluster in self.__clusters.itervalues() if len(cluster) >= min_size)

    def export_bad_clusters(self,
                            error_rate=None,
                            min_size=None,
                            pv_threshold=0.01,
                            out=".",
                            gzip=False):
        if min_size is None:
            min_size = self.__min_size

        if error_rate is None:
            error_rate = self.error_rate()

        ext = ".gz" if gzip else ""

        result = []
        mkdir_p(out)
        for id, cluster in self.__clusters.iteritems():
            pv = cluster.pvalue_upper(error_rate)
            if len(cluster) >= min_size and pv < pv_threshold:
                size = len(cluster)
                merrors = cluster.max_error()
                cluster.export(out + "/bad_cluster__%s__size__%d__merrors__%d__pv_%.5f.fa%s" % (id, size, merrors, pv, ext))
                result.append(id)

        return result

    def plot_estimation_of_max_error_distribution(self,
                                                  error_rate=None,
                                                  out="max_errors",
                                                  pv_threshold=0.01,
                                                  title="",
                                                  annotate=False,
                                                  points=False,
                                                  format=None,
                                                  cis=False):
        import numpy as np
        import matplotlib.pyplot as plt

        if error_rate is None:
            error_rate = self.error_rate()
            if title:
                title += "$\hat{\lambda} = %.3f$" % error_rate
        else:
            if title:
                title += "$\lambda = %.3f$" % error_rate

        f, ax = initialize_plot()

        # TODO take actual cluster size
        read_len = np.median([cluster.length() for cluster in self.__clusters.itervalues()])

        # print "Median read len ", read_len

        def quantile_max_poisson(q, size):
            from scipy.stats import poisson
            import math
            l = read_len
            prob_error = float(error_rate) / float(l) / 3
            lam = float(prob_error) * float(size)
            return poisson.ppf(math.pow(q, 1. / l), mu=lam)

        sizes = [len(cluster) for cluster in self.__clusters.itervalues()]
        max_errors = [cluster.max_error() for cluster in self.__clusters.itervalues()]
        pvalues = [cluster.pvalue_upper(error_rate) for cluster in self.__clusters.itervalues()]
        ids = [cluster.name for cluster in self.__clusters.itervalues()]

        sizes = np.array(sizes)
        max_errors = np.array(max_errors)
        pvalues = np.array(pvalues)
        ids = np.array(ids)
        good = pvalues >= pv_threshold

        if points:
            ax.plot(sizes[good], max_errors[good], "bo", color="blue")
            ax.plot(sizes[~good], max_errors[~good], "bo", color="red")
        else:
            # ax.set_axis_bgcolor('blue')
            f, ax = initialize_plot(figsize=(7.5, 6))
            hexbinblue(ax, sizes, max_errors)

        if annotate:
            for name, size, max_error in zip(ids[~good], sizes[~good], max_errors[~good]):
                plt.annotate(name,
                             xy=(size, max_error),
                             xytext=(size, max_error))

        sizes = range(1, int(ax.get_xlim()[1]) + 1)

        level = 0.05
        q_upper = [quantile_max_poisson(1. - level / 2, size) for size in sizes]
        q50 = [quantile_max_poisson(0.50, size) for size in sizes]
        q_lower = [quantile_max_poisson(level / 2, size) for size in sizes]
        if cis:
            ax.plot(sizes, q_upper, "--", color="black", linewidth=0.5)
            ax.plot(sizes, q_lower, "--", color="black", linewidth=0.5)
            ax.plot(sizes, q50, ":", color="black", linewidth=0.5)

        ax.set_xlabel("Cluster size")
        ax.set_ylabel("maximum secondary votes along profile")

        if title:
            plt.title(title)

        save_plot(out, format=format)

    def plot_cluster_error_profile(self,
                                   out="error_profile",
                                   bins=25,
                                   title="",
                                   min_size=None,
                                   legend=False,
                                   format=None):
        import matplotlib.pyplot as plt
        # import seaborn as sns
        import numpy as np
        from matplotlib import patches

        if min_size is None:
            min_size = self.__min_size

        f, ax = initialize_plot()
        errors01 = self.__errors01()
        values, bins = np.histogram(errors01, bins=bins, range=(0., 1.))

        # HARD fix:
        values[-1:] = np.median(values)
        values[-2:-1] = np.median(values[:-2])

        widths = bins[1:] - bins[:-1]
        width = widths[0]

        width_in_bases = width * 450  # TODO Use len from self.

        values / width_in_bases

        xs = bins[:-1] + widths / 2.

        cdr1_start = 0.25
        cdr1_end = 0.3
        cdr2_start = 0.41
        cdr2_end = 0.54
        cdr3_start = 0.8
        cdr3_end = 0.86

        cdr_mask = ((cdr1_start < xs) & (xs < cdr1_end)) | ((cdr2_start < xs) & (xs < cdr2_end)) | ((cdr3_start < xs) & (xs < cdr3_end))
        cdr_values = values[cdr_mask]
        cdr_bins = bins[:-1][cdr_mask]

        eps = 1. / len(values) / 10
        cdr_color = "#EFBEBE"

        ax.bar(left=cdr_bins + eps,
               height=cdr_values,
               width=width - 2 * eps,
               align="edge",
               # edgecolor='red',
               color=cdr_color,
               label="CDRs")

        plt.xlim(0, 1)
        ax.set_xlabel("Relative position in read")
        ax.set_ylabel("Number of errors")

        if legend:
            handles, labels = ax.get_legend_handles_labels()
            plt.legend(handles, labels)

        ax.bar(left=bins[:-1] + eps,
               height=values,
               width=widths - 2 * eps,
               align="edge",
               # edgecolor='cornflowerblue',
               color='cornflowerblue')

        max_value = ax.get_ylim()[1]
        plt.gca().add_patch(patches.Rectangle((cdr1_start, 0), cdr1_end - cdr1_start, max_value, facecolor=cdr_color, lw=0))
        plt.gca().add_patch(patches.Rectangle((cdr2_start, 0), cdr2_end - cdr2_start, max_value, facecolor=cdr_color, lw=0))
        plt.gca().add_patch(patches.Rectangle((cdr3_start, 0), cdr3_end - cdr3_start, max_value, facecolor=cdr_color, lw=0))

        ax.bar(left=bins[:-1] + eps,
               height=values,
               width=widths - 2 * eps,
               align="edge",
               # edgecolor='cornflowerblue',
               color='cornflowerblue')

        if title:
            plt.title(title)

        save_plot(out, format=format)

    def plot_cluster_error_profile_old(self,
                                       out="error_profile",
                                       bins=25,
                                       title="",
                                       min_size=None,
                                       format=None):
        import matplotlib.pyplot as plt
        # import seaborn as sns
        import numpy as np

        if min_size is None:
            min_size = self.__min_size

        f, ax = initialize_plot()
        errors01 = self.__errors01()
        values, bins = np.histogram(errors01, bins=bins, range=(0., 1.))

        widths = bins[1:] - bins[:-1]
        width = widths[0]

        width_in_bases = width * 450  # TODO Use len from self.

        values / width_in_bases

        xs = bins[:-1] + widths / 2.

        cdr1_start = 0.25
        cdr1_end = 0.3
        cdr2_start = 0.41
        cdr2_end = 0.54
        cdr3_start = 0.8
        cdr3_end = 0.86

        cdr_mask = ((cdr1_start < xs) & (xs < cdr1_end)) | ((cdr2_start < xs) & (xs < cdr2_end)) | ((cdr3_start < xs) & (xs < cdr3_end))
        cdr_values = values[cdr_mask]
        cdr_bins = bins[:-1][cdr_mask]

        eps = 1. / len(values) / 10

        ax.bar(left=cdr_bins + eps,
               height=cdr_values,
               width=width - 2 * eps,
               align="edge",
               # edgecolor='red',
               color='red',
               label="CDRs")

        plt.xlim(0, 1)
        ax.set_xlabel("Relative position in read")
        ax.set_ylabel("Number of errors")

        handles, labels = ax.get_legend_handles_labels()
        plt.legend(handles, labels)

        ax.bar(left=bins[:-1][~cdr_mask] + eps,
               height=values[~cdr_mask],
               width=widths[~cdr_mask] - 2 * eps,
               align="edge",
               # edgecolor='cornflowerblue',
               color='cornflowerblue')

        if title:
            plt.title(title)

        save_plot(out, format=format)

    # TODO make it biplot
    def plot_distribution_of_errors_in_reads(self,
                                             out="ErrorsInReadDistribution",
                                             title="",
                                             max_val=6,
                                             min_size=None,
                                             lam=None,
                                             combine_tail=True,
                                             format=None):
        import numpy as np
        import matplotlib.pyplot as plt
        from scipy.stats import poisson

        nerrors = self.__nerrors_by_read(min_size)

        if lam is None:
            lam = self.error_rate()

        f, ax = initialize_plot()

        bins = np.array(range(max_val + 2)) - 0.5

        labels = ["$%d$" % i for i in xrange(max_val + 1)]
        if combine_tail:
            bins[-1] = np.inf
            labels[-1] = "$\geq %d$" % max_val

        a_heights, _ = np.histogram(nerrors, bins=bins)

        theoretical_distribution = poisson(mu=lam)
        th_probs = theoretical_distribution.cdf(bins[1:]) - theoretical_distribution.cdf(bins[:-1])
        N = len(nerrors)
        b_heights = th_probs * N

        calibrate0 = True
        if calibrate0:
            b_heights = b_heights / b_heights[0] * a_heights[0]

        width = (bins[1] - bins[0]) / 3
        delta_outer = -width / 7
        delta_inner = width / 10

        ax.bar(bins[:-1] + width / 2 + delta_outer,
               a_heights,
               width=width - delta_outer - delta_inner,
               facecolor='cornflowerblue',
               label="Computed frequencies")

        pois_model_name = "Poisson distribution ($\hat{\lambda} = %.2f$)" % lam
        ax.bar(bins[:-1] + width + width / 2 + delta_inner,
               b_heights,
               width=width - delta_outer - delta_inner,
               facecolor='seagreen',
               label=pois_model_name)

        plt.xticks(range(max_val + 1), labels)

        handles, labels = ax.get_legend_handles_labels()
        plt.legend(handles, labels)
        plt.xlim(-width - 2 * max([0, -delta_outer]), max_val + width + 2 * max([0. - delta_outer]))

        if title:
            plt.title(title)

        save_plot(out, format=format)

    def bad_clusters_id(self,
                        error_rate=None,
                        min_size=None,
                        pv_threshold=0.01):
        if error_rate is None:
            error_rate = self.error_rate()

        if min_size is None:
            min_size = self.__min_size

        result = []
        for id, cluster in self.__clusters.iteritems():
            pv = cluster.pvalue_upper(error_rate)
            if len(cluster) >= min_size and pv < pv_threshold:
                result.append(id)

        return result

    def report(self, report, name="constructed_stats"):
        if name not in report.__dict__:
            report.__dict__[name] = {}

        rf = report.__dict__[name]

        rf["error_rate"] = self.error_rate()
        rf["error_rate_estimations"] = self.error_rates().__dict__
        rf["bad_clusters"] = self.bad_clusters_id()


class Report:

    def toText(self, filename):
        with smart_open(filename, "w") as f:
            f.write(str(self))

    def toJson(self, filename,
               indent=4, sort_keys=True,
               **kwargs):
        import json
        with smart_open(filename, "w") as f:
            json.dump(self.__dict__, f,
                      indent=indent, sort_keys=sort_keys,
                      *kwargs)

    def toYaml(self, filename, **kwargs):
        try:
            import yaml
            with smart_open(filename, "w") as f:
                yaml.dump(self.__dict__, f,
                          *kwargs)
        except ImportError:
            print "Module 'yaml' is absent"

    def __str__(self):
        s = ""

        if "reference_based" in self.__dict__:
            rb = self.reference_based
            s += "Reference-based quality measures (with size threshold = %(min_size)d):\n" % rb
            s += "\tSensitivity:\t\t\t\t%(sensitivity)0.4f (%(ref2cons)d / %(reference_size)d)\n" % rb
            s += "\tPrecision:\t\t\t\t%(precision)0.4f (%(cons2ref)d / %(constructed_size)d)\n" % rb
            s += "\tMultiplicity median rate:\t\t%(reference_vs_constructed_size_median_rate)0.4f\n" % rb

            if "jaccard_index" in rb:
                s += "\tClustering similarity measures:\n"
                s += "\t\tJaccard index:\t\t\t%(jaccard_index)0.4f\n" % rb
                s += "\t\tFowlkes-Mallows index:\t\t%(fowlkes_mallows_index)0.4f\n" % rb
                s += "\t\tRand index:\t\t\t%(rand_index)0.4f\n" % rb
                s += "\t\tAdjusted Rand index:\t\t%(adjusted_rand_index)0.4f\n" % rb
                s += "\t\tReference purity:\t\t%(reference_purity)0.4f\n" % rb
                s += "\t\tConstructed purity:\t\t%(constructed_purity)0.4f\n" % rb

                s += "\tOver/under-correction measures:\n"
                s += "\t\tpure constr. clusters:\t%(n_pure_constructed_clusters)d\n" % rb
                s += "\t\tdirty constr. clusters:\t%(n_dirty_constructed_clusters)d\n" % rb
                s += "\t\tconstr. secondary vote reads:\t%(constructed_secondary_votes)d\n" % rb

                s += "\t\tpure ref. clusters:\t%(n_pure_reference_clusters)d\n" % rb
                s += "\t\tdirty ref. clusters:\t%(n_dirty_reference_clusters)d\n" % rb
                s += "\t\tref. secondary vote reads:\t%(reference_secondary_votes)d\n" % rb

                s += "\n"

        if "reference_stats" in self.__dict__:
            st = self.reference_stats
            s += "Reference repertoire statistics:\n"
            s += "\tError rate:\t\t\t\t%(error_rate)0.4f\n" % st
            s += "\n"

        if "constructed_stats" in self.__dict__:
            st = self.constructed_stats
            s += "Constructed repertoire statistics:\n"
            s += "\tError rate:\t\t\t\t%(error_rate)0.4f\n" % st
            s += "\n"

        return s


def json_from_file(filename):
    import json

    with smart_open(filename, "r") as f:
        data = json.load(f)

    return data


def get_plot_various_error_rate_data(dir,
                                     kind="igrec",
                                     what="sensitivity",
                                     woans=False):
    from glob import glob

    suff = "_woans" if woans else ""
    dirnames = glob(dir + "/errate_?.??" + suff)

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
        json = json_from_file(dirname + "/" + kind + "/aimquast/aimquast.json")
        return json["reference_based"][what]

    lambdas = map(extract_lambda, dirnames)
    data = map(extract_data, dirnames)

    lambdas, data = zip(*sorted(zip(lambdas, data)))

    return lambdas, data


def plot_various_error_rate(dir,
                            out="var_error_rate",
                            kind1="igrec", kind2="mixcr", what="sensitivity", woans=False,
                            label1=None, label2=None,
                            title="",
                            format=("png", "pdf", "svg")):
    lambdas, data1 = get_plot_various_error_rate_data(dir, kind=kind1, what=what, woans=woans)
    lambdas, data2 = get_plot_various_error_rate_data(dir, kind=kind2, what=what, woans=woans)
    import matplotlib.pyplot as plt
    import seaborn as sns

    f, ax = initialize_plot()

    if label1 is None:
        label1 = kind1

    if label2 is None:
        label2 = kind2

    sns.set_style("darkgrid")
    plt.plot(lambdas, data1, "b-", color="cornflowerblue", label=label1)
    plt.plot(lambdas, data2, "b-", color="seagreen", label=label2)

    eps = 0.025
    plt.ylim((0. - eps, 1. + eps))

    plt.ylabel(what)
    plt.xlabel("Error rate")

    if title:
        plt.title(title)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc=3)

    save_plot(out, format=format)


def plot_two_rocs(json1, json2,
                  label1="1", label2="2",
                  max_size_threshold=75,
                  out="two_rocs",
                  title="",
                  format=None,
                  plot_second=True):
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns

    json1 = json_from_file(json1)
    json2 = json_from_file(json2)

    sensitivity1 = json1["reference_based"]["__data_sensitivity"]
    sensitivity2 = json2["reference_based"]["__data_sensitivity"]
    precision1 = json1["reference_based"]["__data_precision"]
    precision2 = json2["reference_based"]["__data_precision"]

    if len(sensitivity1) > max_size_threshold:
        sensitivity1 = sensitivity1[:max_size_threshold]
    if len(sensitivity2) > max_size_threshold:
        sensitivity2 = sensitivity2[:max_size_threshold]
    if len(precision1) > max_size_threshold:
        precision1 = precision1[:max_size_threshold]
    if len(precision2) > max_size_threshold:
        precision2 = precision2[:max_size_threshold]

    f, ax = initialize_plot()

    sns.set_style("darkgrid")
    plt.plot(precision1, sensitivity1, "b-", color="cornflowerblue", label=label1)

    if plot_second:
        plt.plot(precision2, sensitivity2, "b-", color="seagreen", label=label2)

    eps = 0.025
    plt.xlim((0 - eps, 1 + eps))
    plt.ylim((0 - eps, 1 + eps))

    plt.ylabel("sensitivity")
    plt.xlabel("precision")

    blue_points = []

    def annotation(size, x, y, color, text=None, xshift=0.01, yshift=-0.04):
        if text is None:
            text = str(size)

        _x = x[size - 1]
        _y = y[size - 1]

        bp = plt.plot(_x, _y, "bo", color=color, label="min size")
        blue_points.append(bp)

        plt.annotate(text, xy=(_x, _y),
                     color=color,
                     xytext=(_x + xshift, _y + yshift))

    for i in [1, 3, 5, 10, 50]:
        annotation(i, precision1, sensitivity1, color="blue")
        if plot_second:
            annotation(i, precision2, sensitivity2, color="green", xshift=-0.04)

    def red_point(x, y):
        reference_trust_cutoff = 5
        if np.isfinite(reference_trust_cutoff):
            plt.plot(x[reference_trust_cutoff - 1], y[reference_trust_cutoff - 1], "bo", color="red", label="Reference size threshold")

    red_point(precision1, sensitivity1)
    if plot_second:
        red_point(precision2, sensitivity2)

    if title:
        plt.title(title)

    handles, labels = ax.get_legend_handles_labels()
    handles = [handles[0], handles[1]]
    labels = [labels[0], labels[1]]

    if not plot_second:
        handles, labels = [handles[0]], [labels[0]]

    ax.legend(handles, labels, loc=3)

    save_plot(out, format=format)
