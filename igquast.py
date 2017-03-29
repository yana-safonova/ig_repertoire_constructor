#!/usr/bin/env python2

import sys
import matplotlib
matplotlib.use('Agg')
import os

current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = current_dir
sys.path.append(igrec_dir + "/py/")
from ash_python_utils import CreateLogger, AttachFileLogger, mkdir_p

sys.path.append(igrec_dir + "/py")

from igquast_impl import Report, reconstruct_rcm, write_rcm, run_consensus_constructor, Repertoire, RepertoireMatch, RcmVsRcm
from igquast_impl import splittering


def parse_command_line():
    import argparse

    def ActionTestFactory(name):
        initial_reads = igrec_dir + "/test_dataset/igquast/%s/input_reads.fa.gz" % name
        import os.path
        if not os.path.isfile(initial_reads):
            return None

        class ActionTest(argparse.Action):

            def __init__(self, option_strings, dest, nargs=None, **kwargs):
                super(ActionTest, self).__init__(option_strings, dest, nargs=0, **kwargs)

            def __call__(self, parser, namespace, values, option_string=None):
                setattr(namespace, "initial_reads", initial_reads)
                setattr(namespace, "output_dir", "igquast_test_%s" % name)
                setattr(namespace, "constructed_repertoire", igrec_dir + "/test_dataset/igquast/%s/igrec/final_repertoire.fa.gz" % name)
                setattr(namespace, "constructed_rcm", igrec_dir + "/test_dataset/igquast/%s/igrec/final_repertoire.rcm" % name)
                setattr(namespace, "reference_repertoire", igrec_dir + "/test_dataset/igquast/%s/repertoire.fa.gz" % name)
                setattr(namespace, "reference_rcm", igrec_dir + "/test_dataset/igquast/%s/repertoire.rcm" % name)

        return ActionTest

    parser = argparse.ArgumentParser(description="IgQUAST: a tool for adaptive immune repertoires quality assessment",
                                     epilog="Report bugs to <igtools_support@googlegroups.com>")

    def add_test(name, key=None, display_name=None):
        if key is None:
            key = "--test-" + name
        if display_name is None:
            display_name = name

        test_action = ActionTestFactory(name)
        if test_action is not None:  # Test dataset exists
            parser.add_argument(key,
                                action=test_action,
                                help="Running on %s dataset" % display_name)

    add_test("test", key="--test")
    add_test("age1")
    add_test("age3")
    add_test("flu")
    add_test("presentation")

    def add_selector(group,
                     key,
                     default,
                     nokey=None,
                     var=None,
                     help=None,
                     help_true=None,
                     help_false=None):
        if var is None:
            var = key.replace("--", "").replace("-", "_")
        if help_true is None:
            help_true = "enable " + help
        if help_false is None:
            help_false = "disable " + help
        if default is True:
            help_true += " (default)"
        else:
            help_false += " (default)"
        if nokey is None:
            nokey = key.replace("--", "--no-")

        selector = group.add_mutually_exclusive_group(required=False)
        selector.add_argument(key,
                              dest=var,
                              action="store_true",
                              help=help_true)
        selector.add_argument(nokey,
                              dest=var,
                              action="store_false",
                              help=help_false)
        group.set_defaults(**{var: default})

        return group

    input = parser.add_argument_group("Input")
    input.add_argument("--initial-reads", "-s",
                       type=str,
                       default="",
                       help="FASTA/FASTQ file with initial reads, empty string for non-providing (default: <empty>)")
    input.add_argument("--constructed-repertoire", "-c",
                       type=str,
                       default="",
                       help="constructed repertoire file")
    input.add_argument("--constructed-rcm", "-C",
                       type=str,
                       default="",
                       help="constructed RCM file, empty string for non-providing (default: <empty>)")
    input.add_argument("--reference-repertoire", "-r",
                       type=str,
                       default="",
                       help="reference reperoire file, empty string for non-providing (default: <empty>)")
    input.add_argument("--reference-rcm", "-R",
                       type=str,
                       default="",
                       help="reference repertoire RCM, empty for non-providing (default: <empty>)")
    add_selector(input,
                 "--reconstruct",
                 default=False,
                 help_true="reconstruct missing repertoire files if possible",
                 help_false="do not reconstruct missing repetoire files")

    output = parser.add_argument_group("Output")
    output.add_argument("--output-dir", "-o",
                        type=str,
                        help="output dir for results")

    output.add_argument("--figure-format", "-F",
                        type=str,
                        default="png,pdf,svg",
                        help="format(s) for produced figures separated by commas, empty for non-producing figures (default: %(default)s)")

    output.add_argument("--json",
                        type=str,
                        help="file for JSON report output (default: <output_dir>/igquast.json)")
    output.add_argument("--text",
                        type=str,
                        help="file for text report output (default: <output_dir>/igquast.txt)")

    scenarios = parser.add_argument_group("Performed scenarios")
    add_selector(scenarios,
                 "--repertoire-to-repertoire-matching",
                 default=True,
                 help_true="perform repertoire-to-repertoire matching",
                 help_false="do not perform repertoire-to-repertoire matching")

    add_selector(scenarios,
                 "--partition-based",
                 default=True,
                 help="partition-based metrics and plots")

    add_selector(scenarios,
                 "--export-bad-clusters",
                 default=False,
                 help_true="export bad clusters during reference-free analysis",
                 help_false="do not export bad clusters during reference-free analysis")

    add_selector(scenarios,
                 "--reference-free",
                 default=False,
                 help="reference-free metrics")

    scenarios.add_argument("--experimental",
                           default=False,
                           action="store_true",
                           help=argparse.SUPPRESS)

    scenarios.add_argument("--page-mode",
                           action="store_true",
                           default=False,
                           help=argparse.SUPPRESS)

    params = parser.add_argument_group("Additional algorithm parameters")
    params.add_argument("--reference-size-cutoff",
                        default=5,
                        help="reference size cutoff (default: %(default)d)")
    params.add_argument("--threads", "-t",
                        type=int,
                        default=16,
                        help="the number of parallel threads used for repertoire-to-reprtoire matching (default: %(default)d)")
    params.add_argument("--tau",
                        type=int,
                        default=6,
                        # help="maximal distance for repertoire-to-repertoire matching (default: %(default)d)")
                        help=argparse.SUPPRESS)

    args = parser.parse_args()

    args.reference_trash_cutoff = args.reference_trust_cutoff = args.reference_size_cutoff
    assert 0 < args.reference_trash_cutoff <= args.reference_trust_cutoff

    if args.output_dir is None:
        parser.print_help()
        sys.exit(1)

    args.log = args.output_dir + "/igquast.log"

    if args.json is None:
        args.json = args.output_dir + "/igquast.json"

    if args.text is None:
        args.text = args.output_dir + "/igquast.txt"

    if args.export_bad_clusters:
        args.reference_free = True

    args.rcm_based = args.reference_free or args.partition_based

    args.reference_free_dir = args.output_dir + "/reference_free"
    args.reference_based_dir = args.output_dir + "/reference_based"

    args.figure_format = [fmt.strip() for fmt in args.figure_format.strip().split(",")]
    args.figure_format = [fmt for fmt in args.figure_format if fmt in ["svg", "pdf", "png"]]

    return args


def main(args):
    if args.reconstruct:
        if args.initial_reads and args.constructed_repertoire and not args.constructed_rcm:
            log.info("Try to reconstruct repertoire RCM file...")
            rcm = reconstruct_rcm(args.initial_reads, args.constructed_repertoire, threads=args.threads)
            args.constructed_rcm = args.output_dir + "/constructed.rcm"
            write_rcm(rcm, args.constructed_rcm)

        if args.initial_reads and not args.constructed_repertoire and args.constructed_rcm:
            log.info("Try to reconstruct repertoire sequence file...")
            args.constructed_repertoire = args.output_dir + "/constructed.fa.gz"
            run_consensus_constructor(rcm_file=args.constructed_rcm,
                                      initial_reads=args.initial_reads,
                                      output_file=args.constructed_repertoire,
                                      threads=args.threads)

        if args.initial_reads and args.reference_repertoire and not args.reference_rcm and args.rcm_based:
            log.info("Try to reconstruct reference RCM file...")
            rcm = reconstruct_rcm(args.initial_reads, args.reference_repertoire, threads=args.threads)
            args.reference_rcm = args.output_dir + "/reference.rcm"
            write_rcm(rcm, args.reference_rcm)

        if args.initial_reads and not args.reference_repertoire and args.reference_rcm and args.rcm_based:
            log.info("Try to reconstruct reference repertoire sequence file...")
            args.reference_repertoire = args.output_dir + "/reference.fa.gz"
            run_consensus_constructor(rcm_file=args.reference_rcm,
                                      initial_reads=args.initial_reads,
                                      output_file=args.reference_repertoire,
                                      threads=args.threads)

    if args.initial_reads and args.constructed_repertoire and args.constructed_rcm and args.rcm_based:
        rep = Repertoire(args.constructed_rcm, args.initial_reads, args.constructed_repertoire)

    report = Report()
    report.min_size = args.reference_size_cutoff

    def ref_free_plots(rep, name, dir):
        if args.figure_format:
            mkdir_p(dir)

            rep.plot_distribution_of_errors_in_reads(out=dir + "/%s_distribution_of_errors_in_reads" % name,
                                                     format=args.figure_format)
            rep.plot_estimation_of_max_error_distribution(out=dir + "/%s_max_error_scatter" % name,
                                                          format=args.figure_format)
            if args.page_mode:
                ymax = 0.041
            else:
                ymax = 0

            for i in range(5):
                cluster = rep.largest(i)
                cluster.plot_profile(out=dir + "/%s_cluster_discordance_profile_largest_%d" % (name, i + 1),
                                     format=args.figure_format,
                                     ymax=ymax)
                cluster.plot_profile(out=dir + "/%s_cluster_error_profile_largest_%d" % (name, i + 1),
                                     discordance=False,
                                     format=args.figure_format,
                                     ymax=ymax)

            rep.plot_profile(out=dir + "/%s_cluster_discordance_profile" % name,
                             format=args.figure_format,
                             ymax=ymax)
            rep.plot_profile(out=dir + "/%s_cluster_error_profile" % name,
                             discordance=False,
                             format=args.figure_format,
                             ymax=ymax)

        if args.export_bad_clusters:
            mkdir_p(dir)
            rep.export_bad_clusters(out=dir + "/bad_%s_clusters/" % name)
        rep.report(report, "%s_stats" % name)

    if args.initial_reads and args.constructed_repertoire and args.constructed_rcm and args.reference_free:
        ref_free_plots(rep, "constructed", args.reference_free_dir)

    if args.initial_reads and args.reference_repertoire and args.reference_rcm and args.reference_free:
        rep_ideal = Repertoire(args.reference_rcm, args.initial_reads, args.reference_repertoire)
        ref_free_plots(rep_ideal, "reference", args.reference_free_dir)

    if args.constructed_repertoire and args.reference_repertoire:
        res = RepertoireMatch(args.constructed_repertoire,
                              args.reference_repertoire,
                              tmp_file=None,
                              max_tau=args.tau,
                              reference_trash_cutoff=args.reference_trash_cutoff,
                              reference_trust_cutoff=args.reference_trust_cutoff,
                              log=log,
                              threads=args.threads)

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

                res.plot_octoplot(out=args.reference_based_dir + "/distance_distribution",
                                  format=args.figure_format)

            res.plot_min_cluster_size_choose(out=args.reference_based_dir + "/sensitivity_precision",
                                             format=args.figure_format)

            res.plot_error_pos_dist(out=args.reference_based_dir + "/error_position_distribution",
                                    format=args.figure_format)

            res.plot_reference_vs_constructed_size(out=args.reference_based_dir + "/cluster_abundances_scatterplot",
                                                   format=args.figure_format, marginals=False)

            if args.experimental:
                res.plot_reference_vs_constructed_size(out=args.reference_based_dir + "/cluster_abundances_scatterplot_hexes",
                                                       points=False,
                                                       format=args.figure_format, marginals=False)

            res.plot_multiplicity_distributions(out=args.reference_based_dir + "/abundance_distributions_log",
                                                ylog=True,
                                                format=args.figure_format)

            res.plot_multiplicity_distributions(out=args.reference_based_dir + "/abundance_distributions",
                                                ylog=False,
                                                format=args.figure_format)

    if args.constructed_rcm and args.reference_rcm and args.partition_based:
        rcm2rcm = RcmVsRcm(args.constructed_rcm,
                           args.reference_rcm)

        rcm2rcm.report(report, "rcm_stats_all_clusters")

        size = args.reference_size_cutoff

        rcm2rcm_large = rcm2rcm.prune_copy(size, 1)

        rcm2rcm_large.report(report)

        if args.figure_format:
            mkdir_p(args.reference_based_dir)
            if args.page_mode:
                ymax = 17000
            else:
                ymax = 0

            rcm2rcm.plot_purity_distribution(out=args.reference_based_dir + "/constructed_purity_distribution", format=args.figure_format, ymax=ymax)
            rcm2rcm.plot_discordance_distribution(out=args.reference_based_dir + "/constructed_discordance_distribution", format=args.figure_format, ymax=ymax)
            if args.experimental:
                rcm2rcm.plot_purity_distribution(out=args.reference_based_dir + "/constructed_purity_distribution_ylog", format=args.figure_format, ylog=True, ymax=ymax)
                rcm2rcm.plot_discordance_distribution(out=args.reference_based_dir + "/constructed_discordance_distribution_ylog", format=args.figure_format, ylog=True, ymax=ymax)

            rcm2rcm.plot_purity_distribution(out=args.reference_based_dir + "/reference_purity_distribution", format=args.figure_format, constructed=False, ymax=ymax)
            rcm2rcm.plot_discordance_distribution(out=args.reference_based_dir + "/reference_discordance_distribution", format=args.figure_format, constructed=False, ymax=ymax)
            if args.experimental:
                rcm2rcm.plot_purity_distribution(out=args.reference_based_dir + "/reference_purity_distribution_ylog", format=args.figure_format, constructed=False, ylog=True, ymax=ymax)
                rcm2rcm.plot_discordance_distribution(out=args.reference_based_dir + "/reference_discordance_distribution_ylog", format=args.figure_format, constructed=False, ylog=True, ymax=ymax)

            rcm2rcm_large.plot_purity_distribution(out=args.reference_based_dir + "/constructed_purity_distribution_large", format=args.figure_format, ymax=ymax)
            rcm2rcm_large.plot_discordance_distribution(out=args.reference_based_dir + "/constructed_discordance_distribution_large", format=args.figure_format, ymax=ymax)
            if args.experimental:
                rcm2rcm_large.plot_purity_distribution(out=args.reference_based_dir + "/constructed_purity_distribution_large_ylog", format=args.figure_format, ylog=True, ymax=ymax)
                rcm2rcm_large.plot_discordance_distribution(out=args.reference_based_dir + "/constructed_discordance_distribution_large_ylog", format=args.figure_format, ylog=True, ymax=ymax)

            rcm2rcm_large.plot_purity_distribution(out=args.reference_based_dir + "/reference_purity_distribution_large", format=args.figure_format, constructed=False, ymax=ymax)
            rcm2rcm_large.plot_discordance_distribution(out=args.reference_based_dir + "/reference_discordance_distribution_large", format=args.figure_format, constructed=False, ymax=ymax)
            if args.experimental:
                rcm2rcm_large.plot_discordance_distribution(out=args.reference_based_dir + "/reference_discordance_distribution_large_ylog", format=args.figure_format, constructed=False, ylog=True, ymax=ymax)
                rcm2rcm_large.plot_purity_distribution(out=args.reference_based_dir + "/reference_purity_distribution_large_ylog", format=args.figure_format, constructed=False, ylog=True, ymax=ymax)

            if args.experimental:
                rcm2rcm.plot_majority_secondary(out=args.reference_based_dir + "/constructed_majority_secondary", format=args.figure_format)
                rcm2rcm.plot_size_nomajority(out=args.reference_based_dir + "/constructed_size_nomajority", format=args.figure_format)
                rcm2rcm.plot_majority_secondary(out=args.reference_based_dir + "/reference_majority_secondary", format=args.figure_format, constructed=False)
                rcm2rcm.plot_size_nomajority(out=args.reference_based_dir + "/reference_size_nomajority", format=args.figure_format, constructed=False)
                rcm2rcm_large.plot_majority_secondary(out=args.reference_based_dir + "/constructed_majority_secondary_large", format=args.figure_format)
                rcm2rcm_large.plot_size_nomajority(out=args.reference_based_dir + "/constructed_size_nomajority_large", format=args.figure_format)
                rcm2rcm_large.plot_majority_secondary(out=args.reference_based_dir + "/reference_majority_secondary_large", format=args.figure_format, constructed=False)
                rcm2rcm_large.plot_size_nomajority(out=args.reference_based_dir + "/reference_size_nomajority_large", format=args.figure_format, constructed=False)

    if args.constructed_rcm and args.reference_rcm and args.constructed_repertoire and args.reference_repertoire and args.experimental:
        splittering(rcm2rcm, rep, args, report)

    log.info(report)

    if args.text:
        report.toText(args.text)

    if args.json:
        report.toJson(args.json)


def SupportInfo(log):
    log.info("\nIn case you have troubles running IgQUAST, "
             "you can write to igtools_support@googlegroups.com.")
    log.info("Please provide us igquast.log file from the output directory.")

if __name__ == "__main__":
    args = parse_command_line()
    mkdir_p(args.output_dir)
    log = CreateLogger("IgQUAST")
    if args.log:
        AttachFileLogger(log, args.log)

    try:
        log.info("Command line: %s" % " ".join(sys.argv))
        main(args)
        log.info("\nThank you for using IgQUAST!")
    except (KeyboardInterrupt):
        log.info("\nIgQUAST has been interrupted!")
    except Exception:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            log.info("\nERROR: Exception caught.")
            SupportInfo(log)
            sys.exit(exc_value)
    except BaseException:
        exc_type, exc_value, _ = sys.exc_info()
        if exc_type == SystemExit:
            sys.exit(exc_value)
        else:
            log.exception(exc_value)
            log.info("\nERROR: Exception caught.")
            SupportInfo(log)
            sys.exit(exc_value)

    log.info("Log was written to " + args.log)
