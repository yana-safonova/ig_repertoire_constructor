from Bio import SeqIO
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# TODO fix fonts
# FROM http://stackoverflow.com/questions/11367736/matplotlib-consistent-font-using-latex


current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = current_dir + "/../"
sys.path.append(igrec_dir + "/py/utils")
sys.path.append(igrec_dir + "/py/pipeline/")
import support
sys.path.append(igrec_dir + "/py/")
from ash_python_utils import idFormatByFileName, smart_open, mkdir_p, FakeLog, memoize, memoize_invalidate

sys.path.append(igrec_dir + "/py")
from ig_compress_equal_clusters import parse_cluster_mult


path_to_igrec = igrec_dir


def safe_log(x):
    from math import log
    try:
        result = log(x)
    except:
        result = float("-inf")
    return result


def auc(sensitivity, precision):
    assert len(precision) == len(sensitivity)

    sensitivity = [sensitivity[0]] + sensitivity + [0]
    precision = [0] + precision + [precision[-1]]

    area = 0
    for xp, xn, yp, yn in zip(precision[:-1], precision[1:], sensitivity[:-1], sensitivity[1:]):
        area += (xn * yp - xp * yn)

    area /= 2.0

    return area


def __opt_fnc(sensitivity, precision, fnc):
    assert len(precision) == len(sensitivity)

    i = max(xrange(len(sensitivity)), key=lambda i: fnc([sensitivity[i] + precision[i]]))
    size = i + 1

    return size, sensitivity[i], precision[i], fnc([sensitivity[i], precision[i]])


def opt_sum(sensitivity, precision):
    return __opt_fnc(sensitivity, precision, sum)


def opt_min(sensitivity, precision):
    return __opt_fnc(sensitivity, precision, min)


def hexbinblue(ax, x, y, gridsize=25, xlim=None, ylim=None):
    import matplotlib.colors as colors
    import matplotlib.pyplot as plt
    plt.hexbin(x, y,
               gridsize=gridsize,
               # bins='log',
               norm=colors.LogNorm(),
               mincnt=1,
               edgecolor="darkblue",
               linewidths=1,
               cmap="Blues")
    plt.axis([xlim[0], xlim[1], ylim[0], ylim[1]])
    plt.colorbar()
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
                    tau=4, k=10, strategy=3, threads=16):
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
            "strategy": strategy,
            "threads": threads}

    support.sys_call("%(swg_cmd)s -r %(reference_file)s \
                     -i %(constructed_file)s -o %(output_file)s \
                     -k %(k)d --tau %(tau)d -A --strategy=%(strategy)s \
                     --threads=%(threads)d" % args,
                     log=log)


def run_consensus_constructor(initial_reads, rcm_file, output_file, log=None, threads=16):
    if log is None:
        log = FakeLog()

    args = {"path": path_to_igrec,
            "cc_cmd": path_to_igrec + '/build/release/bin/ig_component_splitter',
            "initial_reads": initial_reads,
            "output_file": output_file,
            "rcm_file": rcm_file,
            "threads": threads}

    support.sys_call("%(cc_cmd)s -i %(initial_reads)s -R %(rcm_file)s -o %(output_file)s \
                     --threads=%(threads)%d" % args,
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
        if i < len(self.reference_cluster_sizes_unique):
            return self.median_rates_unique[i]
        else:
            print "Ad-hoc fix warning"
            return 0.

    def mean_rate(self, size=1):
        from bisect import bisect_left
        if len(self.reference_cluster_sizes_unique) == 0:
            return 0.

        i = bisect_left(self.reference_cluster_sizes_unique, size)  # the leftmost elt >= size
        if i < len(self.reference_cluster_sizes_unique):
            return self.mean_rates_unique[i]
        else:
            print "Ad-hoc fix warning"
            return 0.

    def plot_reference_vs_constructed_size(self, out, title="", format=None,
                                           points=True,
                                           marginals=False):
        import seaborn as sns

        if points:
            f, ax = initialize_plot()
        else:
            f, ax = initialize_plot(figsize=(7.5, 6))

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

        is_overcorrected = self.constructed_cluster_sizes > self.reference_cluster_sizes

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
                plt.plot(self.reference_cluster_sizes[~is_overcorrected],
                         self.constructed_cluster_sizes[~is_overcorrected], "bo",
                         color="blue",
                         label="clusters", alpha=0.5,
                         markersize=6)
                plt.plot(self.reference_cluster_sizes[is_overcorrected],
                         self.constructed_cluster_sizes[is_overcorrected], "bo",
                         color="red",
                         label="overestimated clusters", alpha=0.5,
                         markersize=6)

            else:
                # ax.set_axis_bgcolor('darkgrey')
                hexbinblue(ax, self.reference_cluster_sizes, self.constructed_cluster_sizes,
                           xlim=(0, uplimit), ylim=(0, uplimit))
                # plt.axes().set_aspect('equal', 'datalim')
                # plt.axes().set_aspect('equal')

            ax.set_xlabel("Reference cluster size")
            ax.set_ylabel("Constructed cluster size")

        ax.plot([0, uplimit], [0, uplimit], "--", linewidth=0.5, color="black", label="$y = x$")
        # ax.plot([0, uplimit], [0, uplimit * self.median_rate(5)])

        # ax.plot(self.reference_cluster_sizes_unique, self.reference_cluster_sizes_unique * self.mean_rates_unique)
        ax.plot(self.reference_cluster_sizes_unique, self.reference_cluster_sizes_unique * self.median_rates_unique,
                label="median")

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc=2)

        if title:
            plt.title(title)

        save_plot(out, format=format)


class Reperoire2RepertoireMatching:

    def __init__(self,
                 constructed_repertoire, reference_repertoire,
                 tmp_file=None, max_tau=4, log=None, threads=16):
        if tmp_file is None:
            import tempfile
            tmp_file = tempfile.mkstemp(suffix=".graph", prefix="igquast_")[1]

        run_ig_matcher2(reference_repertoire, constructed_repertoire,
                        tau=max_tau,
                        output_file=tmp_file, log=log, prefix="",
                        threads=threads)
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
                    tmp_file=None, max_tau=4, log=None, threads=16):
        r2c = Reperoire2RepertoireMatching(constructed_repertoire=constructed_repertoire,
                                           reference_repertoire=reference_repertoire,
                                           tmp_file=tmp_file, max_tau=max_tau, log=log,
                                           threads=threads)
        c2r = Reperoire2RepertoireMatching(constructed_repertoire=reference_repertoire,
                                           reference_repertoire=constructed_repertoire,
                                           tmp_file=tmp_file, max_tau=max_tau, log=log,
                                           threads=threads)

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
                # for i in matches:
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
            bottom = np.ones_like(constructed_h)
            constructed_h -= bottom
            reference_h -= bottom
        else:
            bottom = None

        ax.bar(bins[:-1] - width + delta_outer,
               reference_h,
               width=width - delta_outer - delta_inner,
               facecolor='cornflowerblue',
               label="Reference abundances")

        ax.bar(bins[:-1] + delta_inner,
               constructed_h,
               width=width - delta_outer - delta_inner,
               facecolor='seagreen',
               label="Constructed abundances")

        # plt.xticks(range(max_val + 1), labels)
        xlim = ax.get_xlim()
        xlim = (-width, xlim[1])
        ylim = ax.get_ylim()
        ylim = (0, ylim[1])
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel("Cluszer size")
        ax.set_ylabel("#clusters")

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
    def __read_fa_cluster_ids(filename):
        with smart_open(filename) as f:
            result = map(lambda rec: parse_cluster_mult(rec.id)[0], SeqIO.parse(f, idFormatByFileName(filename)))

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

    @memoize
    def extra_clusters(self, cons_size=5, ref_size=5):
        assert len(self.constructed_ids) == len(self.rep2rep.constructed_abundances) == len(self.constructed_new[:, 0])
        return [id for id, ab, match in zip(self.constructed_ids, self.rep2rep.constructed_abundances, self.constructed_new[:, 0]) if ab >= cons_size and match < ref_size]

    @memoize
    def extra_clusters_ref_sizes(self, cons_size=5, ref_size=5):
        return [match for ab, match in zip(self.rep2rep.constructed_abundances, self.constructed_new[:, 0]) if ab >= cons_size and match < ref_size]

    @memoize
    def missed_clusters(self, cons_size=5, ref_size=5):
        assert len(self.reference_ids) == len(self.rep2rep.reference_abundances) == len(self.reference_new[:, 0])
        return [id for id, ab, match in zip(self.reference_ids, self.rep2rep.reference_abundances, self.reference_new[:, 0]) if ab >= cons_size and match < ref_size]

    @memoize
    def missed_clusters_cons_sizes(self, cons_size=5, ref_size=5):
        return [match for ab, match in zip(self.rep2rep.reference_abundances, self.reference_new[:, 0]) if ab >= cons_size and match < ref_size]

    def __init__(self,
                 constructed_repertoire, reference_repertoire,
                 tmp_file=None, max_tau=10,
                 reference_trash_cutoff=-float("inf"),
                 reference_trust_cutoff=float("inf"),
                 log=None,
                 threads=16):
        import numpy as np

        self.reference_ids = self.__read_fa_cluster_ids(reference_repertoire)
        self.constructed_ids = self.__read_fa_cluster_ids(constructed_repertoire)

        self.reference_repertoire = self.__read_fa(reference_repertoire)
        self.constructed_repertoire = self.__read_fa(constructed_repertoire)

        self.rep2rep = Reperoire2RepertoireMatching.bidirection(constructed_repertoire=constructed_repertoire,
                                                                reference_repertoire=reference_repertoire,
                                                                max_tau=max_tau,
                                                                tmp_file=tmp_file,
                                                                log=log,
                                                                threads=threads)

        self.rep2rep.check(log=log,
                           reference_repertoire=self.reference_repertoire,
                           constructed_repertoire=self.constructed_repertoire)

        assert reference_trash_cutoff <= reference_trust_cutoff

        ref_len = len(self.rep2rep.reference_abundances)
        cons_len = len(self.rep2rep.constructed_abundances)

        reference = np.full((ref_len, max_tau + 1), 0, dtype=int)
        constructed = np.full((cons_len, max_tau + 1), 0, dtype=int)

        reference_new = np.full((ref_len, max_tau + 1), 0, dtype=int)
        constructed_new = np.full((cons_len, max_tau + 1), 0, dtype=int)

        ref_edge = np.full((ref_len, max_tau + 1), -1, dtype=int)
        cons_edge = np.full((cons_len, max_tau + 1), -1, dtype=int)

        reference_sum = np.full((ref_len, max_tau + 1), 0, dtype=int)
        constructed_sum = np.full((cons_len, max_tau + 1), 0, dtype=int)

        assert len(self.rep2rep.constructed2reference) == len(self.constructed_repertoire)
        assert len(self.rep2rep.reference2constructed) == len(self.reference_repertoire)

        for i, pairs in enumerate(self.rep2rep.constructed2reference):
            for j, d in pairs:
                # assert self.__hamming(self.constructed_repertoire[i], self.reference_repertoire[j]) == d

                reference_abundance = self.rep2rep.reference_abundances[j]
                constructed_abundance = self.rep2rep.constructed_abundances[i]

                reference_sum[j, d] = reference_sum[j, d] + constructed_abundance
                constructed_sum[i, d] = constructed_sum[i, d] + reference_abundance

                reference_new[j, d] = max(reference_new[j, d], constructed_abundance)
                constructed_new[i, d] = max(constructed_new[i, d], reference_abundance)

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
                reference_new[j, d] = max(reference_new[j, d], reference_new[j, d - 1])
                reference_sum[j, d] = reference_sum[j, d] + reference_sum[j, d - 1]

            for i in xrange(cons_len):
                constructed[i, d] = max(constructed[i, d], constructed[i, d - 1])
                constructed_new[i, d] = max(constructed_new[i, d], constructed_new[i, d - 1])
                constructed_sum[i, d] = constructed_sum[i, d] + constructed_sum[i, d - 1]

        self.cons_edge = cons_edge
        self.ref_edge = ref_edge
        self.constructed = constructed
        self.reference = reference
        self.constructed_new = constructed_new
        self.reference_new = reference_new

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

        rb["precision"] = self.precision(size)
        rb["sensitivity"] = self.sensitivity(size)
        rb["cons2ref"] = self.cons2ref(size)
        rb["ref2cons"] = self.ref2cons(size)

        rb["constructed_size"] = self.constructed_size(size)
        rb["reference_size"] = self.reference_size(size)

        rb["reference_vs_constructed_size_median_rate"] = float(self.M2MDATA.median_rate(size))
        rb["reference_vs_constructed_size_mean_rate"] = float(self.M2MDATA.mean_rate(size))  # TODO check type
        rb["supernode_mode_reference_vs_constructed_error_rate_estimation"] = -safe_log(float(self.M2MDATA.median_rate(size)))
        rb["reference_vs_constructed_constructed_clusters_sizes"] = map(int, list(self.M2MDATA.constructed_cluster_sizes))
        rb["reference_vs_constructed_reference_clusters_sizes"] = map(int, list(self.M2MDATA.reference_cluster_sizes))
        cf = 1.0
        rb["reference_vs_constructed_overcorrected"] = sum(self.M2MDATA.constructed_cluster_sizes > cf * self.M2MDATA.reference_cluster_sizes)
        rb["reference_vs_constructed_num"] = len(self.M2MDATA.constructed_cluster_sizes)
        is_large = self.M2MDATA.reference_cluster_sizes >= size
        is_overcorrected = self.M2MDATA.constructed_cluster_sizes > cf * self.M2MDATA.reference_cluster_sizes
        rb["reference_vs_constructed_overcorrected_large"] = sum(self.M2MDATA.constructed_cluster_sizes[is_large] > cf * self.M2MDATA.reference_cluster_sizes[is_large])
        rb["reference_vs_constructed_large_num"] = sum(is_large)
        import numpy as np
        rb["reference_vs_constructed_overcorrected_large_median_rate"] = np.median(self.M2MDATA.constructed_cluster_sizes[is_large & is_overcorrected] / self.M2MDATA.reference_cluster_sizes[is_large & is_overcorrected])
        is_large50 = self.M2MDATA.reference_cluster_sizes >= 50
        rb["reference_vs_constructed_overcorrected_large50"] = sum(self.M2MDATA.constructed_cluster_sizes[is_large50] > cf * self.M2MDATA.reference_cluster_sizes[is_large50])
        rb["reference_vs_constructed_large50_num"] = sum(is_large50)
        rb["reference_vs_constructed_overcorrected_large_median_rate50"] = np.median(self.M2MDATA.constructed_cluster_sizes[is_large50 & is_overcorrected] / self.M2MDATA.reference_cluster_sizes[is_large50 & is_overcorrected])
        rb["reference_vs_constructed_overcorrected_large_mean_rate50"] = np.mean(self.M2MDATA.constructed_cluster_sizes[is_large50 & is_overcorrected] / self.M2MDATA.reference_cluster_sizes[is_large50 & is_overcorrected])
        rb["reference_vs_constructed_overcorrected_large_mean_rate"] = np.mean(self.M2MDATA.constructed_cluster_sizes[is_large & is_overcorrected] / self.M2MDATA.reference_cluster_sizes[is_large & is_overcorrected])

        precision, sizes = self.__get_data(what="precision")
        sensitivity, _ = self.__get_data(what="sensitivity")
        precision, sensitivity, sizes = list(precision), list(sensitivity), list(sizes)
        rb["__data_sizes"] = sizes
        rb["__data_sensitivity"] = sensitivity
        rb["__data_precision"] = precision

        rb["AUC"] = auc(sensitivity, precision)
        rb["opt_min_size"], rb["opt_min_sensitivity"], rb["opt_min_precision"], rb["opt_min"] = opt_min(sensitivity, precision)
        rb["opt_sum_size"], rb["opt_sum_sensitivity"], rb["opt_sum_precision"], rb["opt_sum"] = opt_sum(sensitivity, precision)

        ref2cons5, taus, _ = self.__get_measure_for_plotting(size=5, what="ref2cons")
        cons2ref5, _, _ = self.__get_measure_for_plotting(size=5, what="cons2ref")
        rb["__data_taus"] = list(taus)
        rb["__data_ref2cons5"] = list(ref2cons5)
        rb["__data_cons2ref5"] = list(cons2ref5)

        rb["error_positions_reversed"] = self.__error_positions(reversed=True)
        rb["error_positions"] = self.__error_positions(reversed=False)

        rb["missed_clusters"] = self.missed_clusters()
        rb["extra_clusters"] = self.extra_clusters()
        rb["missed_clusters_sizes"] = self.missed_clusters_cons_sizes()
        rb["extra_clusters_ref_sizes"] = self.extra_clusters_ref_sizes()

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

        f, ax = initialize_plot(font_scale=1.5)

        if max_tau is None:
            max_tau = self.max_tau

        data = self.__get_measure_for_plotting(size=size, what=what, differential=differential, max_tau=max_tau)
        measures = np.array(data[0])
        taus = np.array(data[1])
        labels = data[2]

        width = 0.9

        if what in ["sensitivity", "ref2cons"]:
            color = "cornflowerblue"
        elif what in ["precision", "cons2ref"]:
            color = "orangered"

        ax.bar(taus + 0.5 - width / 2, measures, width=width,
               facecolor=color,
               label="Actual frequencies")
        plt.xticks(taus + 0.5,
                   labels)

        plt.xlim(0, max(taus) + 1)
        plt.ylim(0, sum(measures))

        if what in ["sensitivity", "ref2cons"]:
            # plt.title("Distribution of distance from reference to constructed sequences")
            plt.title("Distance from reference to constructed")
        elif what in ["precision", "cons2ref"]:
            # plt.title("Distribution of distance from constructed to reference sequences")
            plt.title("Distance from constructed to reference")

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
            ax.set_ylabel("# of errors")
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
                                     legend=False,
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

        if legend:
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
                if cluster.strip() == "":
                    cluster = None
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

        self.ids = ids

    @memoize_invalidate
    def fix_nones(self):
        def fix_nones(v, prefix):
            for i in xrange(len(v)):
                if v[i] is None:
                    v[i] = prefix + str(i)

        fix_nones(self.clustering1, "__none__clustering1__")
        fix_nones(self.clustering2, "__none__clustering2__")

    @memoize
    def fixed_nones(self):
        from copy import deepcopy
        result = deepcopy(self)
        result.fix_nones()
        return result

    @memoize
    def indices(self):
        fixed = self.fixed_nones()

        indices = clustering_similarity_indices(fixed.clustering1, fixed.clustering2)
        indices.constructed_purity = purity(fixed.clustering1, fixed.clustering2)
        indices.reference_purity = purity(fixed.clustering2, fixed.clustering1)

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
        # return votes(self.clustering1, self.clustering2) if constructed else votes(self.clustering2, self.clustering1)
        return self.votes_dict(constructed).values()

    @memoize
    def votes_dict(self, constructed=True):
        return votes_dict(self.clustering1, self.clustering2) if constructed else votes_dict(self.clustering2, self.clustering1)

    @memoize
    def purity(self, constructed=True):
        import numpy as np

        votes = self.votes(constructed)
        majority_votes = [vote[0] for vote in votes]
        # secondary_votes = [vote[1] for vote in votes]
        sizes = [sum(vote) for vote in votes]

        return np.array(majority_votes, dtype=float) / np.array(sizes, dtype=float)

    @memoize
    def discordance(self, constructed=True):
        import numpy as np

        votes = self.votes(constructed)
        # majority_votes = [vote[0] for vote in votes]
        secondary_votes = [vote[1] for vote in votes]
        sizes = [sum(vote) for vote in votes]

        return np.array(secondary_votes, dtype=float) / np.array(sizes, dtype=float)

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
            ax.set_xlabel("Primary votes")
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
            ax.set_ylabel("Cluster size - first votes")

            save_plot(out, format=format)
        except BaseException as ex:
            print ex

    def plot_purity_distribution(self, out, **kwargs):
        self.plot_purity_discordance_distribution(out, what="purity", **kwargs)

    def plot_discordance_distribution(self, out, **kwargs):
        self.plot_purity_discordance_distribution(out, what="purity", **kwargs)

    def plot_purity_discordance_distribution(self, out,
                                             format=None, constructed=True, ylog=False, what=None,
                                             xmax=None, ymax=0):
        import seaborn as sns

        assert what in ["purity", "discordance"]

        f, ax = initialize_plot()

        if xmax is None:
            xmax = 1 if what == "purity" else 0.5


        data = self.discordance(constructed) if what == "discordance" else self.purity(constructed)
        try:
            sns.distplot(data, kde=False, bins=25, ax=ax)
            ax.set_xlabel("Discordance" if what == "discordance" else "Purity")
            ax.set_ylabel("#clusters")
            ax.set_xlim((0, xmax))
            if ylog:
                plt.yscale("log", nonposy="clip")
            else:
                ymax = max(ymax, len(data))
                ax.set_ylim((0, ymax))

            save_plot(out, format=format)
        except BaseException as ex:
            print ex


def reconstruct_rcm(initial_reads, repertoire,
                    tmp_file_matcher=None, tmp_file_reads=None,
                    taus=(1, 2, 4, 8, 12),
                    fallback_to_exhaustive_mode=False,
                    log=None,
                    threads=16):
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
                                                       log=log, max_tau=tau,
                                                       threads=threads)

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


def votes_dict(X, Y):
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

    votes = {cluster: sorted(cluster_content.itervalues(), reverse=True) for cluster, cluster_content in cluster.iteritems()}

    for vote in votes.itervalues():
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

    # nuc2idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    # reads = [str(read.seq) for read in reads]

    # In [5]: np.fromstring("ACGT", count=4, dtype=np.uint8)
    # Out[5]: array([65, 67, 71, 84], dtype=uint8)

    mx = np.zeros(shape=(l, 4), dtype=int)
    for read in reads:
        read = np.fromstring(str(read.seq), dtype=np.uint8, count=l)
        # assert len(read) == l
        mx[read == 65, 0] += 1
        mx[read == 67, 1] += 1
        mx[read == 71, 2] += 1
        mx[read == 84, 3] += 1
        # assert sum(read == 65) + sum(read == 67) + sum(read == 71) + sum(read == 84) == len(read)

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
    def nbases(self):
        return sum([len(read) for read in self.__reads])

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
        ms = majority_secondary(self.__reads)
        l = self.length()
        return ms[0][:l], ms[1][:l]

    @memoize
    def errors_by_positions(self):
        return len(self) - self.majority_secondary()[0]

    @memoize
    def secondary_votes_by_positions(self):
        return self.majority_secondary()[1]

    @memoize
    def nerrors_by_read(self):
        import numpy as np

        center = self.centroid()
        center = np.fromstring(str(center.seq), dtype=np.uint8, count=len(center.seq))
        errors_by_read = []

        for read in self.__reads:
            read = np.fromstring(str(read.seq), dtype=np.uint8, count=len(read.seq))
            l = min(len(read), len(center))
            error_by_read = np.sum(read[:l] != center[:l])
            errors_by_read.append(error_by_read)

        return errors_by_read

    @memoize
    def errors01(self):
        result = []
        errors_by_positions = self.errors_by_positions()

        l = self.length()
        for i, nerrors in enumerate(errors_by_positions):
            if nerrors:
                result += [float(i) / float(l)] * nerrors

        return result

    @memoize
    def max_error(self):
        return max(self.errors_by_positions()) if self.length() else 0

    @memoize
    def max_secondary_vote(self):
        return max(self.secondary_votes_by_positions()) if self.length() else 0

    @memoize
    def discordance(self):
        return self.max_secondary_vote(self) / float(len(self))

    @memoize
    def max_ppf(self, q, size, error_rate):
        from scipy.stats import poisson
        import math

        l = self.length()

        prob_error = float(error_rate) / float(l)
        lam = float(prob_error) * float(size)
        return poisson.ppf(math.pow(q, 1. / l), mu=lam)

    @memoize
    def max_cdf(self, x, size, error_rate):
        from scipy.stats import poisson

        l = self.length()

        prob_error = float(error_rate) / float(l)
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
        return sum(self.errors_by_positions())


def join_list_of_lists(a):
    import itertools
    return list(itertools.chain.from_iterable(a))


def sum_npv_list(a):
    import numpy as np

    res = np.zeros(0)
    for v in a:
        if len(v) > len(res):
            res.resize(len(v))
        res[:len(v)] += np.array(v)

    return res


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
            if id in id2read and cluster is not None:  # TODO fix it according to issue #99
                clusters[cluster].append(id2read[id])
                clusters[cluster].name = cluster

        for read in reads:
            cluster, mult = parse_cluster_mult(str(read.description))
            if cluster in clusters:
                clusters[cluster].mult = mult
                clusters[cluster].center = read
                clusters[cluster].name = cluster

        self.clusters = clusters

    @memoize_invalidate
    def __invalidate(self):
        pass

    def largest(self, index=0):
        from collections import defaultdict
        from copy import copy

        result = copy(self)
        len_name = [(len(cluster), cluster.name) for cluster in self.clusters.itervalues()]
        len_name.sort(reverse=True)
        index = min(len(self.clusters), index)
        name = len_name[index][1]
        largest_cluster = self.clusters[name]
        result.clusters = defaultdict(Cluster)
        result.clusters[name] = largest_cluster

        result.__invalidate()

        return result

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
            res.first_len = -safe_log(float(errordist[0]) / N)
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

        return join_list_of_lists(cluster.nerrors_by_read() for cluster in self.clusters.itervalues() if len(cluster) >= min_size)

    @memoize
    def __errors01(self, min_size=None):
        if min_size is None:
            min_size = self.__min_size

        return join_list_of_lists(cluster.errors01() for cluster in self.clusters.itervalues() if len(cluster) >= min_size)

    @memoize
    def __secondary_votes_by_positions(self, min_size=None):
        if min_size is None:
            min_size = self.__min_size

        return sum_npv_list(cluster.secondary_votes_by_positions() for cluster in self.clusters.itervalues() if len(cluster) >= min_size)

    @memoize
    def __errors_by_positions(self, min_size=None):
        if min_size is None:
            min_size = self.__min_size

        return sum_npv_list(cluster.errors_by_positions() for cluster in self.clusters.itervalues() if len(cluster) >= min_size)

    def export_bad_clusters(self,
                            error_rate=None,
                            min_size=None,
                            pv_threshold=0.01,
                            out=".",
                            gzip=True):
        if min_size is None:
            min_size = self.__min_size

        if error_rate is None:
            error_rate = self.error_rate()

        ext = ".gz" if gzip else ""

        result = []
        mkdir_p(out)
        for id, cluster in self.clusters.iteritems():
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
                                                  points=True,
                                                  format=None,
                                                  cis=True):
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
        read_len = np.median([cluster.length() for cluster in self.clusters.itervalues()])

        # print "Median read len ", read_len

        def quantile_max_poisson(q, size):
            from scipy.stats import poisson
            import math
            l = read_len
            prob_error = float(error_rate) / float(l)
            lam = float(prob_error) * float(size)
            return poisson.ppf(math.pow(q, 1. / l), mu=lam)

        sizes = [len(cluster) for cluster in self.clusters.itervalues()]
        max_errors = [cluster.max_error() for cluster in self.clusters.itervalues()]
        pvalues = [cluster.pvalue_upper(error_rate) for cluster in self.clusters.itervalues()]
        ids = [cluster.name for cluster in self.clusters.itervalues()]

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
        ax.set_ylabel("Maximum #errors along read")

        if title:
            plt.title(title)

        save_plot(out, format=format)

    def plot_profile(self,
                     out="discordance_profile",
                     title="",
                     discordance=True,
                     min_size=None,
                     legend=False,
                     format=None,
                     ymax=0):
        import matplotlib.pyplot as plt
        # import seaborn as sns
        import numpy as np
        from matplotlib import patches

        if min_size is None:
            min_size = self.__min_size

        f, ax = initialize_plot()

        if discordance:
            values = self.__secondary_votes_by_positions()
        else:
            values = self.__errors_by_positions()

        l = len(values)
        bins = np.arange(l + 1)

        widths = bins[1:] - bins[:-1]
        width = widths[0]

        sumlen = sum([len(cluster) for cluster in self.clusters.itervalues()])

        values = values / float(sumlen)
        # print "Peak", max(values)

        xs = bins[:-1] + widths / 2.

        cdr1_start = 0.25 * l
        cdr1_end = 0.3 * l
        cdr2_start = 0.41 * l
        cdr2_end = 0.54 * l
        cdr3_start = 0.8 * l
        cdr3_end = 0.86 * l

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

        plt.xlim(0, l)
        ax.set_xlabel("Position in constructed cluster")
        if discordance:
            ax.set_ylabel("Discordance")
        else:
            ax.set_ylabel("Error rate")

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
        max_value = max(max_value, ymax)
        plt.gca().add_patch(patches.Rectangle((cdr1_start, 0), cdr1_end - cdr1_start, max_value, facecolor=cdr_color, lw=0))
        plt.gca().add_patch(patches.Rectangle((cdr2_start, 0), cdr2_end - cdr2_start, max_value, facecolor=cdr_color, lw=0))
        plt.gca().add_patch(patches.Rectangle((cdr3_start, 0), cdr3_end - cdr3_start, max_value, facecolor=cdr_color, lw=0))

        ax.bar(left=bins[:-1] + eps,
               height=values,
               width=widths - 2 * eps,
               align="edge",
               # edgecolor='cornflowerblue',
               color='cornflowerblue')
        plt.ylim((0, max_value))

        if title:
            plt.title(title)

        save_plot(out, format=format)

    def plot_distribution_of_errors_in_reads(self,
                                             out="ErrorsInReadDistribution",
                                             title="",
                                             max_val=6,
                                             min_size=None,
                                             lam=None,
                                             combine_tail=True,
                                             additional_space_for_legend=True,
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
        plt.legend(handles, labels, loc="upper right")
        plt.xlim(-width - 2 * max([0, -delta_outer]), max_val + width + 2 * max([0. - delta_outer]))

        if additional_space_for_legend:
            ylim = plt.ylim()
            plt.ylim((ylim[0], 1.2*ylim[1]))

        if title:
            plt.title(title)

        save_plot(out, format=format)

    def bad_clusters_ids(self,
                         error_rate=None,
                         min_size=None,
                         pv_threshold=0.01):
        if error_rate is None:
            error_rate = self.error_rate()

        if min_size is None:
            min_size = self.__min_size

        result = []
        for id, cluster in self.clusters.iteritems():
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
        rf["bad_clusters"] = self.bad_clusters_ids()


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

        min_size = self.min_size
        if "reference_based" in self.__dict__:
            rb = self.reference_based

            if "sensitivity" in rb:
                s += "Reference-based quality measures (with size threshold = %d):\n" % min_size
                s += "\tSensitivity:\t\t\t\t%(sensitivity)0.4f (%(ref2cons)d / %(reference_size)d)\n" % rb
                s += "\tPrecision:\t\t\t\t%(precision)0.4f (%(cons2ref)d / %(constructed_size)d)\n" % rb
                s += "\tAbundances median rate:\t\t\t%(reference_vs_constructed_size_median_rate)0.4f\n" % rb
                s += "\tArea under curve:\t\t\t%(AUC)0.4f\n" % rb
                s += "\tMaximal S + P:\t\t\t\t%(opt_sum_sensitivity)0.4f + %(opt_sum_precision)0.4f = %(opt_sum)0.4f\n" % rb
                s += "\tMaxizing S + P constructed min size:\t%(opt_sum_size)d\n" % rb

                extra_clusters_ref_sizes = rb["extra_clusters_ref_sizes"]
                for size in reversed(range(min_size)):
                    n_ec = sum([1 for sz in extra_clusters_ref_sizes if sz == size])
                    s += "\tExtra clusters with size == %d:\t\t%d\n" % (size, n_ec)

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
                s += "\t\timpure constr. clusters:\t%(n_dirty_constructed_clusters)d\n" % rb
                s += "\t\tconstr. discordance reads:\t%(constructed_secondary_votes)d\n" % rb

                s += "\t\tpure ref. clusters:\t%(n_pure_reference_clusters)d\n" % rb
                s += "\t\timpure ref. clusters:\t%(n_dirty_reference_clusters)d\n" % rb
                s += "\t\tref. discordance reads:\t%(reference_secondary_votes)d\n" % rb

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


def plot_logit(y, X, colors, colormap=False):
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_style("darkgrid")

    # print X.shape
    # print X
    assert X.shape[0] == len(y)
    # print X[:, 1]
    # print X.min()
    # print X[:, 1].min()
    # print X[:, 1].max()
    # clf = LogisticRegression().fit(np.hstack((X, np.ones((X.shape[0], 1)))), y)

    clf = LogisticRegression().fit(X, y)
    score = clf.score(X, y)
    coef = clf.coef_
    intercept = clf.intercept_

    # xx, yy = np.mgrid[X[:, 0].min():X[:, 0].max():.01,
    #                   X[:, 1].min():X[:, 1].max():.01]
    xx, yy = np.meshgrid(np.arange(X[:, 0].min(), X[:, 0].max(), 0.5),
                         np.arange(X[:, 1].min(), X[:, 1].max(), 0.5),
                         sparse=False, indexing='xy')

    grid = np.c_[xx.ravel(), yy.ravel()]
    # grid = np.hstack((grid, np.ones((grid.shape[0], 1))))

    probs = clf.predict_proba(grid)[:, 1].reshape(xx.shape)

    if colormap:
        f, ax = initialize_plot()
        contour = ax.contourf(xx, yy, probs, 25, cmap="RdBu",
                              vmin=0, vmax=1)
        ax_c = f.colorbar(contour)
        ax_c.set_label("$P(y = 1)$")
        ax_c.set_ticks([0, .25, .5, .75, 1])

        ax.scatter(X[:, 0], X[:, 1], c=colors,
                   s=30,
                   cmap="RdBu", vmin=-.2, vmax=1.2,
                   edgecolor="white", linewidth=1)
        plt.xlabel("Cluster size")
        plt.ylabel("Secondary votes")
    else:
        f, ax = initialize_plot()
        # ax.contour(xx, yy, probs, levels=[.5], cmap="Greys", vmin=0, vmax=.6)
        # ax.contour(xx, yy, probs, levels=[.5], cmap="Greys")

        ax.scatter(X[:, 0], X[:, 1], s=20,
                   c=colors,
                   # cmap="RdBu", vmin=-.2, vmax=1.2, alpha=0.7)
                   cmap="RdBu", alpha=0.7)
        plt.xlabel("Cluster size")
        plt.ylabel("Secondary votes")
        xlim = plt.xlim()
        ylim = plt.ylim()
        ax.contour(xx, yy, probs, levels=[.5], cmap="Greys", vmin=0, vmax=.6)
        plt.xlim(0, xlim[1])
        plt.ylim(0, ylim[1])

    return coef, intercept, score


def splittering(rcm2rcm, rep, args, report):
    mp = rcm2rcm.fixed_nones().votes_dict(constructed=True)
    for cluster in rep.clusters.itervalues():
        assert cluster.name in mp
        votes = mp[cluster.name]
        purity = float(votes[0]) / sum(votes)
        cluster.purity = purity

    purities = [cluster.purity for cluster in rep.clusters.itervalues()]
    second_votes = [cluster.max_secondary_vote() for cluster in rep.clusters.itervalues()]
    sizes = [len(cluster) for cluster in rep.clusters.itervalues()]

    import seaborn as sns
    import numpy as np
    import matplotlib.pyplot as plt

    x = np.array(second_votes, dtype=np.float) / np.array(sizes)
    y = np.array(purities)

    sizes = np.array(sizes)

    def F(size):
        ax = sns.regplot(x=x[sizes >= size],
                         y=y[sizes >= size],
                         line_kws={"alpha": 0.7})
        ax.set_xlabel("Secondary votes")
        ax.set_ylabel("Purity")
        eps = 0.05
        plt.ylim(0 - eps, 1 + eps)
        plt.xlim(0 - eps, 0.5 + eps)
        plt.savefig(args.reference_based_dir + "/purity_vs_secondvote_%d.png" % size)
        plt.close()
    map(F, [1, 5, 10, 15, 50])

    os.system(igrec_dir + "/build/release/bin/ig_component_splitter \
              -i %s -o %s -R %s -M %s -V 1 \
              --allow-unassigned=true \
              --recursive=false --flu=false" % (args.initial_reads,
                                                args.output_dir + "/splitted.fa.gz",
                                                args.constructed_rcm,
                                                args.output_dir + "/splitted.rcm"))

    def read_fa_cluster_ids(filename):
        with smart_open(filename) as f:
            result = {parse_cluster_mult(rec.id)[0]: rec for rec in SeqIO.parse(f, idFormatByFileName(filename))}

        return result

    splitted_clusters = read_fa_cluster_ids(args.output_dir + "/splitted.fa.gz")
    constructed_clusters = read_fa_cluster_ids(args.constructed_repertoire)
    repertoire = read_fa_cluster_ids(args.reference_repertoire)
    repertoire = set(str(rec.seq) for rec in repertoire.itervalues() if parse_cluster_mult(rec.description)[1] >= 5)

    score_diffs = []
    cluster_sizes = []
    max_second_votes = []

    for id in constructed_clusters:
        constructed_seq = constructed_clusters[id]

        splitted_seqs = []
        for splitted_id in [id + "X0", id + "X1"]:
            if splitted_id in splitted_clusters:
                splitted_seqs.append(splitted_clusters[splitted_id])

        cluster_size = len(rep.clusters[id])
        initial_score = 0
        assert cluster_size == parse_cluster_mult(constructed_seq.description)[1]
        if cluster_size < 5:
            continue

        constructed_seq = str(constructed_seq.seq)
        initial_score += 1 if constructed_seq in repertoire else -1

        splitted_score = 0
        for seq in splitted_seqs:
            if parse_cluster_mult(seq.description)[1] < 5:
                continue
            seq = str(seq.seq)
            splitted_score += 1 if seq in repertoire else -1

        max_secondary_vote = rep.clusters[id].max_secondary_vote()
        # print initial_score, splitted_score, max_secondary_vote, cluster_size, max_secondary_vote / cluster_size
        score_diffs.append(splitted_score - initial_score)
        cluster_sizes.append(cluster_size)
        max_second_votes.append(max_secondary_vote)

    if "splitting" not in report.__dict__:
        report.splitting = {}

    rs = report.splitting

    rs["score_diffs"] = score_diffs
    rs["cluster_sizes"] = cluster_sizes
    rs["max_second_votes"] = max_second_votes

    score_diffs = np.array(score_diffs)
    cluster_sizes = np.array(cluster_sizes)
    max_second_votes = np.array(max_second_votes)

    def def_color(score):
        if score < 0:
            return "red"
        elif score == 0:
            return "black"
        else:
            return "blue"

    def def_size(score):
        if score <= 1:
            return 15
        elif score == 2:
            return 20
        else:
            return 30

    colors = map(def_color, score_diffs)
    bullet_sizes = map(def_size, score_diffs)
    plt.scatter(x=cluster_sizes, y=max_second_votes,
                s=bullet_sizes,
                c=colors, alpha=0.7)

    plt.xlabel("Cluster size")
    plt.ylabel("Secondary votes")
    xlim = plt.xlim()
    ylim = plt.ylim()
    plt.xlim(0, xlim[1])
    plt.ylim(0, ylim[1])

    save_plot(args.reference_based_dir + "/splitting_efficiency", format=args.figure_format)

    plt.scatter(x=cluster_sizes, y=max_second_votes/cluster_sizes,
                s=20,
                c=colors, alpha=0.7)
    plt.xlabel("Cluster size")
    plt.ylabel("Secondary votes")
    xlim = plt.xlim()
    ylim = plt.ylim()
    eps = 0.05
    plt.xlim(0, xlim[1])
    plt.ylim(0 - eps/10, 1 + eps/10)

    save_plot(args.reference_based_dir + "/splitting_efficiency01", format=args.figure_format)

    import itertools
    # Test strategies
    for cluster_size_threshold, secondvote_threshold in itertools.product([5, 20, 50, 100, 200], [0.5, 0.2, 0.1, 0.05, 0.01]):
        print "score_diffs[cluster_sizes >= %d & max_second_votes/cluster_size >= %f]" % (cluster_size_threshold, secondvote_threshold)
        s = sum(score_diffs[(cluster_sizes >= cluster_size_threshold) & (max_second_votes/cluster_sizes >= secondvote_threshold)])
        print s
        rs["SavedClusters_%d_%0.3f" % (cluster_size_threshold, secondvote_threshold)] = s

    s = sum(score_diffs[0.634*max_second_votes - 0.0076*cluster_sizes > 6.025])
    rs["SavedClustersLogit"] = s
    print "Saved clusters logit", s

    good = score_diffs > 0
    coef, intercept, score = plot_logit(good, np.vstack((cluster_sizes,
                                                         max_second_votes)).T,
                                        colors=colors)

    red_points = len([1 for color in colors if color == "red"])
    black_points = len([1 for color in colors if color == "black"])
    blue_points = len([1 for color in colors if color == "blue"])

    print "Red: ", red_points
    print "Black: ", black_points
    print "Blue: ", blue_points

    rs["BluePoints"] = blue_points
    rs["BlackPoints"] = black_points
    rs["RedPoints"] = red_points

    print "Coef", coef, intercept
    print "Score", score
    rs["LogitScore"] = score
    rs["LogitCoef"] = [coef[0][0], coef[0][1]]
    rs["LigitIntercept"] = intercept[0]
    save_plot(args.reference_based_dir + "/splitting_efficiency_logit", format=args.figure_format)
    plot_logit(good, np.vstack((cluster_sizes,
                                max_second_votes)).T,
               colors=colors,
               colormap=True)

    save_plot(args.reference_based_dir + "/splitting_efficiency_logit_colormap", format=args.figure_format)
