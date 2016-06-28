#!/usr/bin/env python2


from Bio import SeqIO
import sys
import matplotlib.pyplot as plt
import os

CUTAIL = 0


def memoize_simple(method):
    def new_method(self):
        if "__memoise_cache" not in self.__dict__:
            self.__memoise_cache = {}
        if method not in self.__memoise_cache:
            self.__memoise_cache[method] = method(self)
        return self.__memoise_cache[method]
    return new_method


def memoize(method):
    def new_method(self, *args, **kwargs):
        if "__memoise_cache" not in self.__dict__:
            self.__memoise_cache = {}

        call = (method, args, tuple(sorted(kwargs.iteritems())))

        if call not in self.__memoise_cache:
            self.__memoise_cache[call] = method(self, *args, **kwargs)

        return self.__memoise_cache[call]
    return new_method


def memoize_invalidate(method):
    def new_method(self, *args, **kwargs):
        self.__memoise_cache = {}
        return method(self, *args, **kwargs)
    return new_method


# class Test:
#
#     @memoize
#     def test(self):
#         print "test() called"
#         return 1
#
#     @memoize
#     def test_a_plus_b(self, a, b):
#         print "test_a_plus_b() called with", a, b
#         return a + b
#
# test = Test()
# print test.test()
# print test.test()
# print test.test()
# print test.test()
# print test.test_a_plus_b(1, 2)
# print test.test_a_plus_b(1, 2)
# print test.test_a_plus_b(1, 2)
# print test.test_a_plus_b(1, 3)
# print test.test_a_plus_b(1, b=3)
# print test.test_a_plus_b(a=1, b=3)
# print test.test_a_plus_b(a=1, b=3)
# print test.test_a_plus_b(a=1, b=3)
# print test.test_a_plus_b(b=3, a=1)
# sys.exit(0)


# matplotlib.rcParams['mathtext.fontset'] = 'custom'
# matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
# matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
# matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

# matplotlib.rcParams['mathtext.fontset'] = 'stix'
# matplotlib.rcParams['font.family'] = 'STIXGeneral'
# TODO fix fonts
# FROM http://stackoverflow.com/questions/11367736/matplotlib-consistent-font-using-latex


current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = current_dir
sys.path.append(igrec_dir + "/src/ig_tools/python_utils")
sys.path.append(igrec_dir + "/src/python_pipeline/")
import support
sys.path.append(igrec_dir + "/src/extra/ash_python_utils/")
from ash_python_utils import CreateLogger, AttachFileLogger, idFormatByFileName, smart_open, mkdir_p, fastx2fastx

sys.path.append(igrec_dir + "/py")
from ig_compress_equal_clusters import parse_cluster_mult


path_to_ig_simulator = igrec_dir + "/../ig_simulator/"
path_to_mixcr = current_dir + "/src/extra/aimquast/mixcr/"
path_to_igrec = igrec_dir


def parse_final_repertoire_id(id):
    import re
    id = id.strip()

    m = re.match(r"^antibody_(\d+)_multiplicity_(\d+)_copy_(\d+)$", id)

    if m:
        g = m.groups()
        return int(g[0]), int(g[1]), int(g[2])
    else:
        return 1


assert parse_final_repertoire_id("antibody_1_multiplicity_1_copy_1") == (1, 1, 1)


def simulated_repertoire_to_rcm(input_file, rcm_file):
    with open(input_file) as fh, open(rcm_file, "w") as fout:
        for record in SeqIO.parse(fh, "fasta"):
            id = record.description
            cluster = str(parse_final_repertoire_id(id)[0])
            fout.write("%s\t%s\n" % (id, cluster))


def RC(l):
    import random

    S = set(list("ACTG"))
    s = S.difference([l])
    return random.choice(list(s))


def jit_fa_file(input_file, output_file, error_rate=2, random_errors=True,
                min_error=0, erroneous_site_len=10005000, seed=None):
    import numpy as np
    from Bio import Seq
    import random

    output_format = idFormatByFileName(output_file)
    random.seed(seed)

    with smart_open(input_file) as fh, smart_open(output_file, "w") as fout:
        for record in SeqIO.parse(fh, "fasta"):
            n_errors = np.random.poisson(error_rate, 1)[0] if random_errors else error_rate
            if n_errors < min_error:
                n_errors = min_error

            positions = random.sample(range(min(len(record.seq), erroneous_site_len)), n_errors)
            s = list(str(record.seq))
            for pos in positions:
                s[pos] = RC(s[pos])
            record.letter_annotations = {}
            record.seq = Seq.Seq("".join(s))

            if output_format == "fastq":
                record.letter_annotations["phred_quality"] = [random.randint(30, 50) for _ in xrange(len(record))]  # TODO Check it out
                # record.letter_annotations["phred_quality"] = [50] * len(record)  # TODO Check it out

            SeqIO.write(record, fout, output_format)


def simulated_repertoire_to_final_repertoire(input_file, output_file):
    import random

    output_format = idFormatByFileName(output_file)

    with smart_open(input_file) as fh, smart_open(output_file, "w") as fout:
        for record in SeqIO.parse(fh, "fasta"):
            id = record.description
            cluster, size, copy = parse_final_repertoire_id(id)
            if copy == 1:
                record.id = record.description = "cluster___%d___size___%d" % (cluster, size)
                record.letter_annotations = {}

                if output_format == "fastq":
                    record.letter_annotations["phred_quality"] = [random.randint(30, 50) for _ in xrange(len(record))]  # TODO Check it out
                    # record.letter_annotations["phred_quality"] = [50] * len(record)

                SeqIO.write(record, fout, output_format)


def simulate_data(input_file, output_dir, **kwargs):
    mkdir_p(output_dir)
    simulated_repertoire_to_rcm(input_file, "%s/ideal_final_repertoire.rcm" % output_dir)
    simulated_repertoire_to_final_repertoire(input_file, "%s/ideal_final_repertoire.fa" % output_dir)
    jit_fa_file(input_file, "%s/merged_reads.fq" % output_dir, **kwargs)


def run_ig_simulator(output_dir, log=None,
                     chain="HC", num_bases=100, num_mutated=1000, reprtoire_size=5000):
    if log is None:
        log = FakeLog()

    assert chain in ["HC", "LC"]

    args = {"path": path_to_ig_simulator,
            "output_dir": output_dir,
            "chain": chain,
            "num_bases": num_bases,
            "num_mutated": num_mutated,
            "reprtoire_size": reprtoire_size}

    support.sys_call("%(path)s/ig_simulator.py --chain-type %(chain)s --num-bases %(num_bases)d --num-mutated %(num_mutated)d --repertoire-size %(reprtoire_size)d -o %(output_dir)s --skip-drawing" % args,
                     log=log)


def convert_mixcr_output_to_igrec(input_file, output_file):
    with smart_open(input_file) as fh, smart_open(output_file, "w") as fout:
        # Skip header
        fh.next()

        for i, line in enumerate(fh):
            seq, size = line.strip().split()
            size = int(size)
            fout.write(">cluster___%d___size___%d\n" % (i, size))
            fout.write(seq + "\n")


def run_igrec(input_file, output_dir, log=None, tau=4, loci="IGH", additional_args=""):
    if log is None:
        log = FakeLog()

    args = {"path": path_to_igrec,
            "tau": tau,
            "loci": loci,
            "input_file": input_file,
            "output_dir": output_dir,
            "additional_args": additional_args}
    support.sys_call("%(path)s/igrec.py --tau=%(tau)d -t 36 --loci %(loci)s -s %(input_file)s -o %(output_dir)s %(additional_args)s" % args,
                     log=log)


def run_mixcr(input_file, output_dir, log=None, loci="IGH"):
    if log is None:
        log = FakeLog()

    mkdir_p(output_dir)

    if idFormatByFileName(input_file) == "fasta":
        input_file_fq = "%s/input_reads.fq" % output_dir
        fastx2fastx(input_file, input_file_fq)
        input_file = input_file_fq

    args = {"path": path_to_mixcr,
            "compress_eq_clusters_cmd": path_to_igrec + "/py/ig_compress_equal_clusters.py",
            "mixcr_cmd": "java -jar %s/mixcr.jar" % path_to_mixcr,
            "input_file": input_file,
            "output_dir": output_dir,
            "loci": loci}

    support.sys_call("%(mixcr_cmd)s  align -f -g -r %(output_dir)s/align_report.txt --loci %(loci)s --noMerge -OvParameters.geneFeatureToAlign=VTranscript %(input_file)s %(output_dir)s/mixcr.vdjca" % args,
                     log=log)
    support.sys_call("%(mixcr_cmd)s assemble -f -r %(output_dir)s/assemble_report.txt -OassemblingFeatures=\"{FR1Begin:FR4Begin}\" %(output_dir)s/mixcr.vdjca %(output_dir)s/mixcr.clns" % args,
                     log=log)
    support.sys_call("%(mixcr_cmd)s exportClones -pf %(path)s/preset.pf -f --no-spaces %(output_dir)s/mixcr.clns %(output_dir)s/mixcr.txt" % args,
                     log=log)
    convert_mixcr_output_to_igrec("%(output_dir)s/mixcr.txt" % args, "%(output_dir)s/mixcr_uncompressed.fa" % args)
    support.sys_call("%(compress_eq_clusters_cmd)s %(output_dir)s/mixcr_uncompressed.fa %(output_dir)s/mixcr_final.fa" % args)


class FakeLog:

    def info(self, msg):
        print msg

    def warn(self, msg):
        print msg


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


class MultToMultData:

    def __init__(self, constructed_abundances, constructed_sum, tau=0, reversed=False):
        import numpy as np

        constructed2reference = zip(constructed_abundances, constructed_sum[:, tau]) if not reversed else zip(constructed_sum[:, tau], constructed_abundances)
        constructed2reference = [x for x in constructed2reference if min(x[0], x[1]) > 0]
        constructed2reference.sort(key=lambda x: x[1])
        constructed_cluster_sizes, reference_cluster_sizes = map(lambda x: np.array(x, dtype=float), zip(*constructed2reference))

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

        i = bisect_left(self.reference_cluster_sizes_unique, size)  # the leftmost elt >= size
        assert i < len(self.reference_cluster_sizes_unique)
        return self.median_rates_unique[i]

    def mean_rate(self, size=1):
        from bisect import bisect_left

        i = bisect_left(self.reference_cluster_sizes_unique, size)  # the leftmost elt >= size
        assert i < len(self.reference_cluster_sizes_unique)
        return self.mean_rates_unique[i]

    def plot_reference_vs_constructed_size(self, out, title="", format=None):
        import seaborn as sns

        f, ax = initialize_plot()

        uplimit = max(self.reference_cluster_sizes)

        def round_up(number, ndigits=0):
            from math import ceil
            p = 10 ** ndigits
            return ceil(number * p) / p
        uplimit = round_up(uplimit, -2)

        g = sns.JointGrid(x=self.reference_cluster_sizes,
                          y=self.constructed_cluster_sizes,
                          xlim=(0, uplimit), ylim=(0, uplimit),
                          ratio=5).set_axis_labels("Reference cluster size", "Constructed cluster size")

        g = g.plot_joint(sns.plt.scatter)

        g.plot_marginals(sns.distplot, kde=True, color=".5")

        ax = g.ax_joint
        ax.plot([0, uplimit], [0, uplimit], "--", linewidth=0.5, color="black", label="$y = x$")
        # ax.plot([0, uplimit], [0, uplimit * self.median_rate(5)])

        # ax.plot(self.reference_cluster_sizes_unique, self.reference_cluster_sizes_unique * self.mean_rates_unique)
        ax.plot(self.reference_cluster_sizes_unique, self.reference_cluster_sizes_unique * self.median_rates_unique,
                label="rate median smoothing")

        g.ax_joint.collections[0].set_label("clusters")
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

    def check(self, log=None):
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
                ref.append((j, matches))

        return cons, ref

    def plot_multiplicity_distributions(self,
                                        out,
                                        bins=25,
                                        ylog=True,
                                        format=format):
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

    def __init__(self,
                 constructed_repertoire, reference_repertoire,
                 tmp_file=None, max_tau=4,
                 reference_trash_cutoff=-float("inf"),
                 reference_trust_cutoff=float("inf"),
                 log=None):
        import numpy as np

        self.rep2rep = Reperoire2RepertoireMatching.bidirection(constructed_repertoire=constructed_repertoire,
                                                                reference_repertoire=reference_repertoire,
                                                                max_tau=max_tau,
                                                                log=log)

        self.rep2rep.check(log=log)

        assert reference_trash_cutoff <= reference_trust_cutoff

        ref_len = len(self.rep2rep.reference_abundances)
        cons_len = len(self.rep2rep.constructed_abundances)

        reference = np.full((ref_len, max_tau + 1), 0, dtype=int)
        constructed = np.full((cons_len, max_tau + 1), 0, dtype=int)
        reference_sum = np.full((ref_len, max_tau + 1), 0, dtype=int)
        constructed_sum = np.full((cons_len, max_tau + 1), 0, dtype=int)

        for i, pairs in enumerate(self.rep2rep.constructed2reference):
            for j, d in pairs:
                reference_abundance = self.rep2rep.reference_abundances[j]
                constructed_abundance = self.rep2rep.constructed_abundances[i]

                reference_sum[j, d] = reference_sum[j, d] + constructed_abundance
                constructed_sum[i, d] = constructed_sum[i, d] + reference_abundance

                if reference_abundance >= reference_trust_cutoff:
                    reference_abundance = float("inf")
                if reference_abundance < reference_trash_cutoff:
                    reference_abundance = -float("inf")

                min_abundance = min(constructed_abundance, reference_abundance)
                reference[j, d] = max(reference[j, d], min_abundance)
                constructed[i, d] = max(constructed[i, d], min_abundance)

        for d in xrange(1, max_tau + 1):
            for j in xrange(ref_len):
                reference[j, d] = max(reference[j, d], reference[j, d - 1])
                reference_sum[j, d] = reference_sum[j, d] + reference_sum[j, d - 1]

            for i in xrange(cons_len):
                constructed[i, d] = max(constructed[i, d], constructed[i, d - 1])
                constructed_sum[i, d] = constructed_sum[i, d] + constructed_sum[i, d - 1]

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

    def __get_measure_for_plotting(self,
                                   size=1,
                                   what="sensitivity",
                                   differential=False,
                                   max_tau=4):
        assert what in ["sensitivity", "fdr", "precision", "recall", "ref2cons", "cons2ref"]

        max_tau = min(max_tau, self.max_tau)

        func = lambda tau: getattr(self, what)(size=size, tau=tau)

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
                                   max_tau=4,
                                   format=None):
        import numpy as np
        import matplotlib.pyplot as plt

        f, ax = initialize_plot()

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

        if what in ["sensitivity", "ref2cons"]:
            plt.title("Distribution of distance from reference to constructed sequences")
        elif what in ["precision", "cons2ref"]:
            plt.title("Distribution of distance from reference to constructed sequences")

        ax.set_xlabel("distance")
        ax.set_ylabel("#clusters")
        if what not in ["ref2cons", "cons2ref"]:
            plt.ylim(0, 1)

        save_plot(out, format=format)

    def plot_octoplot(self, out,
                      sizes=(1, 3, 5, 10),
                      differential=True,
                      max_tau=4,
                      format=format):
        import numpy as np
        import matplotlib.pyplot as plt

        f, ax = initialize_plot(figsize=(8, 4))

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
                    facecolor='cornflowerblue')
            plt.xticks(taus + 0.5,
                       labels)
            plt.title("precision, size = %d" % size)

        plt.tight_layout()

        save_plot(out, format=format)

    def plot_min_cluster_size_choose(self,
                                     what_x="precision",
                                     what_y="sensitivity",
                                     max_size_threshold=100,
                                     tau=0,
                                     out="fdr_sensitivity",
                                     title="",
                                     format=format):
        import numpy as np
        import matplotlib.pyplot as plt
        import seaborn as sns

        assert tau <= self.max_tau

        assert what_x in ["sensitivity", "fdr", "precision", "recall"]
        assert what_y in ["sensitivity", "fdr", "precision", "recall"]

        func_x = lambda size: getattr(self, what_x)(size=size, tau=tau)
        func_y = lambda size: getattr(self, what_y)(size=size, tau=tau)

        f, ax = initialize_plot()

        sizes = np.array(range(1, max_size_threshold + 1))

        x = np.array(map(func_x, sizes))
        y = np.array(map(func_y, sizes))

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

            bp = plt.plot(_x, _y, "bo", color="blue", label="min_size")
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

    def __init__(self, rcm1, rcm2, fix_nones=True, size=1):
        rcm1, rcm2 = parse_rcm(rcm1), parse_rcm(rcm2)

        sizes1 = self.rcm_count_cluster_sizes(rcm1)
        sizes2 = self.rcm_count_cluster_sizes(rcm2)

        ids = list(set(rcm1.keys() + rcm2.keys()))
        N = len(ids)

        self.clustering1 = []
        self.clustering2 = []

        for id in ids:
            self.clustering1.append(rcm1[id] if id in rcm1 else None)
            self.clustering2.append(rcm2[id] if id in rcm2 else None)

        for i in xrange(N):
            if sizes1[self.clustering1[i]] < size:
                self.clustering1[i] = None
            if sizes2[self.clustering2[i]] < size:
                self.clustering2[i] = None

        def fix_nones(v, prefix):
            for i in xrange(len(v)):
                if v[i] is None:
                    v[i] = prefix + str(i)

        self.clustering1_none = self.clustering1
        self.clustering2_none = self.clustering2
        if fix_nones:
            fix_nones(self.clustering1, "__none__clustering1__")
            fix_nones(self.clustering2, "__none__clustering2__")

        res = clustering_simularity_indices(self.clustering1, self.clustering2)
        for k, v in res.__dict__.iteritems():
            self.__dict__[k] = v

        self.constructed_purity = purity(self.clustering1, self.clustering2)
        self.reference_purity = purity(self.clustering2, self.clustering1)

    def report(self, report):
        if "reference_based" not in report.__dict__:
            report.reference_based = {}

        rb = report.reference_based

        for k, v in self.__dict__.iteritems():
            rb[k] = v

    @memoize
    def votes(self):
        return votes(self.clustering1_none, self.clustering2_none)

    def plot_majority_secondary(self, out, format):
        import numpy as np
        import matplotlib.pyplot as plt
        import seaborn as sns

        f, ax = initialize_plot()

        votes = self.votes()
        majority_votes = [vote[0] for vote in votes]
        secondary_votes = [vote[1] for vote in votes]
        # sizes = [sum(vote) for vote in votes]

        ax.plot(majority_votes, secondary_votes, "bo", label="clusters")
        ax.set_xlabel("Majority votes")  # Primary
        ax.set_ylabel("Secondary votes")

        save_plot(out, format=format)

    def plot_purity_distribution(self, out, format):
        import numpy as np
        import matplotlib.pyplot as plt
        import seaborn as sns

        f, ax = initialize_plot()

        votes = self.votes()
        majority_votes = [vote[0] for vote in votes]
        # secondary_votes = [vote[1] for vote in votes]
        sizes = [sum(vote) for vote in votes]

        purity = np.array(majority_votes, dtype=float) / np.array(sizes, dtype=float)

        sns.distplot(purity, kde=False, bins=25, ax=ax)
        ax.set_ylabel("Purity")
        ax.set_xlim((0, 1))

        save_plot(out, format=format)


def reconstruct_rcm(initial_reads, repertoire,
                    tmp_file_matcher=None, tmp_file_reads=None,
                    taus=(1, 2, 4, 8, 12, 16, 20, 24),
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


def clustering_simularity_indices(X, Y):
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

    class Struct():
        pass

    result = Struct()
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

    majority_votes = sum(max(cluster_content.itervalues()) for cluster_content in cluster.itervalues())

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


def initialize_plot(figsize=(6, 6)):
    import matplotlib
    import matplotlib.pyplot as plt
    import seaborn as sns
    matplotlib.rcParams.update({'font.size': 16})

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
        mkdir_p(dirname)

    if 'png' in format:
        plt.savefig(plot_name + '.png', format='png')

    if 'pdf' in format:
        from matplotlib.backends.backend_pdf import PdfPages
        pp = PdfPages(plot_name + '.pdf')
        pp.savefig()
        pp.close()

    if 'svg' in format:
        plt.savefig(plot_name + '.svg', format='svg')

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
    def nerrors_by_position(self):
        reads = [str(read.seq) for read in self.__reads]
        from collections import defaultdict
        # TODO get rid off dict here
        center = self.centroid()
        errors_by_position = defaultdict(int)

        for read in reads:
            for i in xrange(min(len(read), len(center) - CUTAIL)):
                if read[i] != center[i]:
                    errors_by_position[i] += 1
        errors_by_position = dict2list(errors_by_position)

        return errors_by_position

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
        nerrors_by_position = self.nerrors_by_position()
        return max(nerrors_by_position) if nerrors_by_position else 0

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

        nerrors = self.__nerrors_by_read()

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
                                                  format=None):
        import numpy as np
        import matplotlib.pyplot as plt

        if error_rate is None:
            error_rate = self.error_rate()
            title += "$\hat{\lambda} = %.3f$" % error_rate
        else:
            title += "$\lambda = %.3f$" % error_rate

        f, ax = initialize_plot()

        # TODO take actual cluster size
        read_len = np.median([cluster.length() for cluster in self.__clusters.itervalues()])

        # print "Median read len ", read_len

        def Q(q, size):
            from scipy.stats import poisson
            import math
            l = read_len
            prob_error = float(error_rate) / float(l)
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

        ax.plot(sizes[good], max_errors[good], "bo", color="blue")
        ax.plot(sizes[~good], max_errors[~good], "bo", color="red")

        if annotate:
            for name, size, max_error in zip(ids[~good], sizes[~good], max_errors[~good]):
                plt.annotate(name,
                             xy=(size, max_error),
                             xytext=(size, max_error))

        sizes = range(1, int(ax.get_xlim()[1]) + 1)

        level = 0.05
        q_upper = [Q(1. - level / 2, size) for size in sizes]
        q50 = [Q(0.50, size) for size in sizes]
        q_lower = [Q(level / 2, size) for size in sizes]
        ax.plot(sizes, q_upper, "--", color="black", linewidth=0.5)
        ax.plot(sizes, q_lower, "--", color="black", linewidth=0.5)
        ax.plot(sizes, q50, ":", color="black", linewidth=0.5)

        ax.set_xlabel("Cluster size")
        ax.set_ylabel("maximum mismatches along profile")

        if title:
            plt.title(title)

        save_plot(out, format=format)

    def plot_cluster_error_profile(self,
                                   out="error_profile",
                                   bins=25,
                                   title="",
                                   min_size=None,
                                   format=format):
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
        cdr_bins = bins[cdr_mask]

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
                                             format=format):
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


def parse_command_line(description="aimQUAST"):
    import argparse

    class ActionTest(argparse.Action):

        def __init__(self, option_strings, dest, nargs=None, **kwargs):
            super(ActionTest, self).__init__(option_strings, dest, nargs=0, **kwargs)

        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, "initial_reads", "aimquast_test_dataset/merged_reads.fq")
            setattr(namespace, "output_dir", "aimquast_test")
            setattr(namespace, "constructed_repertoire", "aimquast_test_dataset/igrec_bad/final_repertoire.fa")
            setattr(namespace, "constructed_rcm", "aimquast_test_dataset/igrec_bad/final_repertoire.rcm")
            setattr(namespace, "reference_repertoire", "aimquast_test_dataset/ideal_final_repertoire.fa")
            setattr(namespace, "reference_rcm", "aimquast_test_dataset/ideal_final_repertoire.rcm")
            setattr(namespace, "json", "aimquast_test/aimquast.json")
            setattr(namespace, "yaml", "aimquast_test/aimquast.yaml")
            setattr(namespace, "text", "aimquast_test/aimquast.txt")

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("--test",
                        action=ActionTest,
                        default="",
                        help="Running of test dataset")
    parser.add_argument("--initial-reads", "-s",
                        type=str,
                        default="",
                        help="FASTA/FASTQ file with initial reads, empty string for non-providing (default: <empty>)")
    parser.add_argument("--constructed-repertoire", "-c",
                        type=str,
                        default="",
                        help="constructed repertoire file")
    parser.add_argument("--constructed-rcm", "-C",
                        type=str,
                        default="",
                        help="constructed RCM file, empty string for non-providing (default: <empty>)")
    parser.add_argument("--reference-repertoire", "-r",
                        type=str,
                        default="",
                        help="reference reperoire file, empty string for non-providing (default: <empty>)")
    parser.add_argument("--reference-rcm", "-R",
                        type=str,
                        default="",
                        help="reference repertoire RCM, empty for non-providing (default: <empty>)")
    parser.add_argument("--output-dir", "-o",
                        type=str,
                        help="output dir for results")
    parser.add_argument("--reference-size-cutoff",
                        default=5,
                        help="reference size cutoff")
    parser.add_argument("--json",
                        help="file for JSON output")
    parser.add_argument("--yaml",
                        help="file for YAML output")
    parser.add_argument("--text",
                        help="file for text output")
    parser.add_argument("--figure-format", "-F",
                        default="png",
                        help="format(s) for producing figures, empty for non-producing (default: %(default)s)")

    args = parser.parse_args()

    args.reference_trash_cutoff = args.reference_trust_cutoff = args.reference_size_cutoff
    assert 0 < args.reference_trash_cutoff <= args.reference_trust_cutoff

    assert args.output_dir
    assert args.constructed_repertoire

    args.log = args.output_dir + "/aimquast.log"

    args.reference_free_dir = args.output_dir + "/reference_free"
    args.reference_based_dir = args.output_dir + "/reference_based"

    args.figure_format = [fmt.strip() for fmt in args.figure_format.strip().split(",")]

    return args


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

            s += "\tClustering simularity measures:\n"
            s += "\t\tJaccard index:\t\t\t%(jaccard_index)0.4f\n" % rb
            s += "\t\tFowlkes-Mallows index:\t\t%(fowlkes_mallows_index)0.4f\n" % rb
            s += "\t\tRand index:\t\t\t%(rand_index)0.4f\n" % rb
            s += "\t\tReference purity:\t\t%(reference_purity)0.4f\n" % rb
            s += "\t\tConstructed purity:\t\t%(constructed_purity)0.4f\n" % rb
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

if __name__ == "__main__":
    args = parse_command_line()
    mkdir_p(args.output_dir)

    report = Report()

    log = CreateLogger("aimQUAST")
    if args.log:
        AttachFileLogger(log, args.log)

    if args.initial_reads and args.constructed_repertoire and not args.constructed_rcm:
        log.info("Try to reconstruct repertoire RCM file...")
        rcm = reconstruct_rcm(args.initial_reads, args.constructed_repertoire)
        args.constructed_rcm = args.output_dir + "/constructed.rcm"
        write_rcm(rcm, args.constructed_rcm)

    if args.initial_reads and args.reference_repertoire and not args.reference_rcm:
        log.info("Try to reconstruct reference RCM file...")
        rcm = reconstruct_rcm(args.initial_reads, args.reference_repertoire)
        args.reference_rcm = args.output_dir + "/reference.rcm"
        write_rcm(rcm, args.reference_rcm)

    if args.initial_reads and args.reference_repertoire and args.reference_rcm:
        mkdir_p(args.reference_free_dir)

        rep_ideal = Repertoire(args.reference_rcm, args.initial_reads, args.reference_repertoire)

        if args.figure_format:
            rep_ideal.plot_cluster_error_profile(out=args.reference_free_dir + "/reference_cluster_error_profile",
                                                 format=args.figure_format)
            rep_ideal.plot_distribution_of_errors_in_reads(out=args.reference_free_dir + "/reference_distribution_of_errors_in_reads",
                                                           format=args.figure_format)
            rep_ideal.plot_estimation_of_max_error_distribution(out=args.reference_free_dir + "/reference_estimation_of_max_error_distribution",
                                                                format=args.figure_format)

        rep_ideal.export_bad_clusters(out=args.reference_free_dir + "/bad_reference_clusters/")
        rep_ideal.report(report, "reference_stats")

    if args.initial_reads and args.constructed_repertoire and args.constructed_rcm:
        mkdir_p(args.reference_free_dir)

        rep = Repertoire(args.constructed_rcm, args.initial_reads, args.constructed_repertoire)

        if args.figure_format:
            rep.plot_cluster_error_profile(out=args.reference_free_dir + "/construced_cluster_error_profile",
                                           format=args.figure_format)
            rep.plot_distribution_of_errors_in_reads(out=args.reference_free_dir + "/constructed_distribution_of_errors_in_reads",
                                                     format=args.figure_format)
            rep.plot_estimation_of_max_error_distribution(out=args.reference_free_dir + "/constructed_estimation_of_max_error_distribution",
                                                          format=args.figure_format)

        rep.export_bad_clusters(out=args.reference_free_dir + "/bad_constructed_clusters/")
        rep.report(report, "constructed_stats")

    if args.constructed_repertoire and args.reference_repertoire:
        res = RepertoireMatch(args.constructed_repertoire,
                              args.reference_repertoire,
                              tmp_file=None,
                              max_tau=4,
                              reference_trash_cutoff=args.reference_trash_cutoff,
                              reference_trust_cutoff=args.reference_trust_cutoff,
                              log=log)

        res.report(report)

        if args.figure_format:
            mkdir_p(args.reference_based_dir)

            for size in [1, 3, 5, 10]:
                res.plot_sensitivity_precision(what="ref2cons",
                                               out=args.reference_based_dir + "/reference_to_constructed_distance_distribution_size_%d" % size,
                                               size=size, differential=True,
                                               format=args.figure_format)

                res.plot_sensitivity_precision(what="cons2ref",
                                               out=args.reference_based_dir + "/constructed_to_reference_distance_distribution_size_%d" % size,
                                               size=size, differential=True,
                                               format=args.figure_format)

                res.plot_octoplot(out=args.reference_based_dir + "/octoplot",
                                  format=args.figure_format)

            res.plot_min_cluster_size_choose(out=args.reference_based_dir + "/min_cluster_size_choose",
                                             format=args.figure_format)

            res.plot_reference_vs_constructed_size(out=args.reference_based_dir + "/reference_vs_constructed_size",
                                                   format=args.figure_format)

            res.plot_multiplicity_distributions(out=args.reference_based_dir + "/multiplicity_distribution",
                                                format=args.figure_format)

    if args.constructed_rcm and args.reference_rcm:
        clustering_scores = RcmVsRcm(args.constructed_rcm,
                                     args.reference_rcm)

        clustering_scores.report(report)
        if args.figure_format:
            mkdir_p(args.reference_based_dir)
            clustering_scores.plot_majority_secondary(out=args.reference_based_dir + "/majority_secondary", format=args.figure_format)
            clustering_scores.plot_purity_distribution(out=args.reference_based_dir + "/purity_distribution", format=args.figure_format)

    log.info(report)

    if args.text:
        report.toText(args.text)

    if args.json:
        report.toJson(args.json)

    if args.yaml:
        report.toYaml(args.yaml)
        # scores = rcm_vs_rcm(args.constructed_rcm,
        #                     args.reference_rcm, size=100000000)

    # Rands scores (with threshold) Is it possible to EFFICIENTLY compute Rand with
    # threshold???
    #
    # Exclude CDR3 stats
# CMD line
# src/extra/aimquast/aimquast.py -i tmp_dir/merged_reads.fq -c tmp_dir/final_repertoire.fa -o oppo -C tmp_dir/final_repertoire.rcm -r tmp_dir/ideal_final_repertoire.fa -R tmp_dir/ideal_final_repertoire.rcm
# src/extra/aimquast/aimquast.py -i tmp_dir/merged_reads.fq -c tmp_dir/BAD/final_repertoire.fa -o oppo -C tmp_dir/BAD/final_repertoire.rcm -r tmp_dir/ideal_final_repertoire.fa -R tmp_dir/ideal_final_repertoire.rcm
# src/extra/aimquast/aimquast.py -s /ssd/simulated/igrec/vj_finder/cleaned_reads.fa -c /ssd/simulated/igrec/final_repertoire.fa  -C /ssd/simulated/igrec/final_repertoire.rcm -o /ssd/simulated/oppo -r /ssd/simulated/ideal_repertoire.clusters.fa -R /ssd/simulated/ideal_repertoire.rcm
