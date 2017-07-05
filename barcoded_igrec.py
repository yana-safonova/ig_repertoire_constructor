#!/usr/bin/env python2

import os
import os.path
import sys

import igrec

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
pipeline_dir = os.path.join(home_directory, "py/pipeline/")
py_src = os.path.join(home_directory, "py")

sys.path.append(py_src)
sys.path.append(pipeline_dir)
import support


def HelpAndReturn(log, parser, exit_code=0):
    log.write = log.info
    parser.print_help(log)
    sys.exit(exit_code)


def EnsureRequiredParametersSet(params, parser, log):
    if not params.output:
        log.error("Please specify the output directory (-o/--output parameter)")
        HelpAndReturn(log, parser)
    if not params.loci:
        log.error("Please specify loci (-l/--loci parameter)")
        HelpAndReturn(log, parser)


def ParseCommandLineParams(log):
    import argparse

    class ActionTest(argparse.Action):
        def __init__(self, option_strings, dest, nargs=None, **kwargs):
            super(ActionTest, self).__init__(option_strings, dest, nargs=0, **kwargs)

        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, "single_reads", os.path.join(home_directory, "test_dataset/barcodedIgReC_test.fasta"))
            setattr(namespace, "output", "barigrec_test")
            setattr(namespace, "loci", "all")
            setattr(namespace, "no_compilation", True)
            setattr(namespace, "ignore_code_changes", True)

    parser = argparse.ArgumentParser(description="IgReC: an algorithm for construction of antibody repertoire from immunosequencing data",
                            epilog="""
    In case you have troubles running IgReC, you can write to igtools_support@googlegroups.com.
    Please provide us with ig_repertoire_constructor.log file from the output directory.
                            """,
                            add_help=False)

    parser.add_argument("--test",
                        action=ActionTest,
                        help="Run on test dataset")

    req_args = parser.add_argument_group("Input")
    input_args = req_args.add_mutually_exclusive_group(required=False)
    input_args.add_argument("-s",
                            dest="single_reads",
                            type=str,
                            default="",  # FIXME This is only for ace's version of python. Locally it works great w/o it
                            help="Single reads in FASTQ format")

    pair_reads = parser.add_argument_group("Paired-end reads")
    pair_reads.add_argument("-1",
                            type=str,
                            dest="left_reads",
                            help="Left paired-end reads in FASTQ format")
    pair_reads.add_argument("-2",
                            type=str,
                            dest="right_reads",
                            help="Right paired-end reads in FASTQ format")

    out_args = parser.add_argument_group("Output")
    out_args.add_argument("-o", "--output",
                          type=str,
                          default="",
                          help="Output directory. Required")

    vj_align_args = parser.add_argument_group("Algorithm arguments")
    vj_align_args.add_argument("-l", "--loci",
                               type=str,
                               dest="loci",
                               # required=True,
                               help="Loci: IGH, IGK, IGL, IG (all BCRs), TRA, TRB, TRG, TRD, TR (all TCRs) or all. Required")

    vj_align_args.add_argument("--organism",
                               type=str,
                               default="human",
                               dest="organism",
                               help="Organism (human and mouse only are supported for this moment) [default: %(default)s]")

    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument("-t", "--threads",
                               type=int,
                               default=16,
                               dest="num_threads",
                               help="Thread number [default: %(default)d]")
    optional_args.add_argument("--umi-graph-tau",
                               type=int,
                               default=1,
                               dest="umi_graph_tau",
                               help="Maximum allowed mismatches between UMIs, used for UMI correction"
                                    "[default: %(default)d]")
    optional_args.add_argument("--igrec-tau",
                               type=int,
                               default=2,
                               dest="igrec_tau",
                               help="Maximum allowed mismatches between UMI clusters"
                                    "[default: %(default)d]")
    optional_args.add_argument("-f", '--min-fillin',
                               type=float,
                               default=0.6,
                               dest="min_fillin",
                               help='Minimum fill-in of dense subgraphs [default: %(default)f]')
    optional_args.add_argument("-k", "--detect-chimeras",
                               dest="detect_chimeras",
                               action="store_true",
                               help="Detect chimeras after clustering, may take significant amount of time, default behavior")
    optional_args.add_argument("-K", "--no-detect-chimeras",
                               dest="detect_chimeras",
                               action="store_false",
                               help="Do not detect chimeras after clustering")
    optional_args.add_argument("--clustering-thr",
                               type=int,
                               default=20,
                               dest="clustering_threshold",
                               help="Threshold distance to unite clusters")

    dev_args = parser.add_argument_group("Developer arguments")
    dev_args.add_argument("-n", "--min-super-read-size",
                               type=int,
                               default=1000000,
                               dest="min_super_read_size",
                               # help="Minimum super read size [default: %(default)d]")
                               help=argparse.SUPPRESS)
    dev_args.add_argument("-p", "--no-compilation",
                               dest="no_compilation",
                               action="store_true",
                               help=argparse.SUPPRESS)
                               # help="Exclude C++ code compilation from the pipeline")
    dev_args.add_argument("-P", "--do-compilation",
                               dest="no_compilation",
                               action="store_false",
                               help=argparse.SUPPRESS)
                               # help="Exclude C++ code compilation from the pipeline")
    dev_args.add_argument("-c", "--ignore-code",
                               dest="ignore_code_changes",
                               action="store_true",
                               help=argparse.SUPPRESS)
                               # help="Ignore code changes when checking stages dependencies")
    dev_args.add_argument("-C", "--no-ignore-code",
                               dest="ignore_code_changes",
                               action="store_false",
                               help=argparse.SUPPRESS)
                               # help="Ignore code changes when checking stages dependencies")
    dev_args.add_argument("--debug-stages",
                               dest="output_intermediate",
                               action="store_true",
                               help=argparse.SUPPRESS)
                               # help="Output repertoire after each step")
    dev_args.add_argument("--umi-cleavage-length",
                               type=int,
                               default=0,
                               dest="umi_cleavage_length",
                               help=argparse.SUPPRESS)
                               # help="Cleave UMIs by the specified length (testing purposes only) [default: %(default)d]")



    params = parser.parse_args()

    # process help
    if len(sys.argv) == 1:
        HelpAndReturn(log, parser)

    EnsureRequiredParametersSet(params, parser, log)

    return parser, params


def CheckParamsCorrectness(parser, params, log):
    if params.single_reads:
        if not os.path.exists(params.single_reads):
            print "File with reads doesn't exist: ", params.single_reads
            HelpAndReturn(log, parser, 1)
    elif not params.left_reads or not params.right_reads or not os.path.exists(params.left_reads) or not os.path.exists(params.right_reads):
        if not params.left_reads or not params.right_reads:
            print "Both left and right reads should be passed. Otherwise use -s option."
        else:
            print "File with reads doesn't exist: ", params.left_reads if not os.path.exists(params.left_reads) else params.right_reads
        HelpAndReturn(log, parser, 1)


class _StagePrepare:
    MAKEFILE = "Makefile"
    BIN_SEPARATOR = " !! "

    def __init__(self):
        pass

    @staticmethod
    def EnsureExists(path):
        if not os.path.exists(path):
            os.makedirs(path)

    @staticmethod
    def _GetUnusedName(dir, default):
        name = default
        import random
        while os.path.exists(os.path.join(dir, name)):
            name = "tmp" + str(random.randint(0, 1e6))
        return name

    @staticmethod
    def _ReplaceVariables(tmp_file, params, log):
        # Files are rather small, so no problem with first reading them
        lines = [line for line in open(tmp_file, 'r')]
        with open(tmp_file, 'w') as file:
            for idx, line in enumerate(lines):
                if line.count(_StagePrepare.BIN_SEPARATOR) > 1:
                    log.error("Too many binary separators '%s' in the line #%d: '%s'" % (_StagePrepare.BIN_SEPARATOR, idx, line))
                    exit(1)
                if _StagePrepare.BIN_SEPARATOR in line:
                    if params.ignore_code_changes:
                        line = line[:line.find(_StagePrepare.BIN_SEPARATOR)] + '\n'
                    else:
                        line = line.replace(_StagePrepare.BIN_SEPARATOR, " ")
                line = line.replace("%RUN_PATH", home_directory)
                if params.single_reads:
                    line = line.replace("%INPUT", os.path.abspath(params.single_reads))
                if params.left_reads and params.right_reads:
                    line = line.replace("%LEFT_INPUT", os.path.abspath(params.left_reads))
                    line = line.replace("%RIGHT_INPUT", os.path.abspath(params.right_reads))
                line = line.replace("%THREADS", str(params.num_threads))
                line = line.replace("%IGREC_TAU", str(params.igrec_tau))
                line = line.replace("%UMI_GRAPH_TAU", str(params.umi_graph_tau))
                line = line.replace("%LOCI", params.loci)
                line = line.replace("%ORGANISM", params.organism)
                line = line.replace("%DETECT_CHIMERAS", str(params.detect_chimeras))
                line = line.replace("%MIN_SUPER_NODE_SIZE", str(params.min_super_read_size))
                line = line.replace("%MIN_FILLIN", str(params.min_fillin))
                line = line.replace("%UMI_CLEAVAGE_LENGTH", str(params.umi_cleavage_length))
                line = line.replace("%CLUSTERING_THRESHOLD", str(params.clustering_threshold))
                line = line.replace("%DEBUG_STAGES", str(params.output_intermediate))
                if '%' in line:
                    log.error("Not all template variables substituted in the makefile, update barcoded_igrec.py script, line #%d: '%s'" % (idx, line))
                    exit(1)
                file.write(line)

    @classmethod
    def Prepare(self, params, stage_template, log, stage_dest=None, makefile_name=MAKEFILE):
        if not stage_dest:
            stage_dest = stage_template
        import shutil
        dest_dir = os.path.join(params.output, stage_dest)
        self.EnsureExists(dest_dir)
        tmp_file_name = self._GetUnusedName(dest_dir, makefile_name)
        tmp_file = os.path.join(dest_dir, tmp_file_name)
        shutil.copyfile(os.path.join(home_directory, "pipeline_makefiles", stage_template, makefile_name), tmp_file)
        self._ReplaceVariables(tmp_file, params, log)
        if tmp_file_name == makefile_name:
            return
        import filecmp
        makefile = os.path.join(dest_dir, makefile_name)
        if filecmp.cmp(tmp_file, makefile):
            os.remove(tmp_file)
        else:
            log.info("Generated new makefile for " + stage_dest)
            os.remove(makefile)
            shutil.move(tmp_file, makefile)


def InitMakeFiles(params, log):
    _StagePrepare.EnsureExists(params.output)
    _StagePrepare.Prepare(params, ".", log, makefile_name="Makefile_vars")
    _StagePrepare.Prepare(params, "no_compilation" if params.no_compilation else "compilation", log, stage_dest="compilation")
    if params.single_reads:
        _StagePrepare.Prepare(params, "vj_finder_input", log, stage_dest="vj_finder")
    else:
        _StagePrepare.Prepare(params, "merged_reads", log)
        _StagePrepare.Prepare(params, "vj_finder", log)
    _StagePrepare.Prepare(params, "umis", log)
    _StagePrepare.Prepare(params, "umi_clustering", log)
    _StagePrepare.Prepare(params, "intermediate_ig_trie_compressor", log)
    _StagePrepare.Prepare(params, "ig_graph_constructor", log)
    _StagePrepare.Prepare(params, "dense_subgraph_finder", log)
    _StagePrepare.Prepare(params, "ig_consensus_finder", log)
    _StagePrepare.Prepare(params, "final_repertoire", log)
    return os.path.join(params.output, "final_repertoire")


def PrintOutputFiles(params, log):
    log.info("\nBarcodedIgReC output:")
    log.info("  * Cleaned Ig-Seq reads were written to %s" % os.path.join(params.output, "vj_finder", "cleaned_reads.fa"))
    log.info("  * Contaminated (not Ig-Seq) reads were written to %s" % os.path.join(params.output, "vj_finder", "filtered_reads.fa"))
    log.info("  * VJ alignment output was written to %s" % os.path.join(params.output, "vj_finder", "v_alignments.fa"))
    log.info("  * Antibody clusters of final repertoire with read multiplicities were written to %s" % os.path.join(params.output, "final_repertoire", "final_repertoire.fa"))
    log.info("  * Antibody clusters of final repertoire with molecule multiplicities were written to %s" % os.path.join(params.output, "final_repertoire", "final_repertoire_umi.fa.gz"))
    log.info("  * Read-cluster map of final repertoire was written to %s" % os.path.join(params.output, "final_repertoire", "final_repertoire.rcm"))


def SupportInfo(log):
    log.info("\nIn case you have troubles running BarcodedIgReC, "
             "you can write to igtools_support@googlegroups.com.")
    log.info("Please provide us with igrec.log file from the output directory.")


def main():
    log = igrec.CreateLogger()
    parser, params = ParseCommandLineParams(log)
    CheckParamsCorrectness(parser, params, log)
    try:
        if not os.path.exists(params.output):
            os.makedirs(params.output)
        igrec.CreateFileLogger(params, log)
        igrec.PrintCommandLine(log)
        final_dir = InitMakeFiles(params, log)
        # We need freshly compiled version to get actual build info
        if not params.no_compilation:
            support.sys_call("make -C " + os.path.join(os.path.dirname(final_dir), "compilation"), log)
        print "===================Build info==================="
        from py import build_info
        build_info.Log(log)
        print "================================================"
        support.sys_call("make -C " + final_dir, log)
        PrintOutputFiles(params, log)
        log.info("\nThank you for using BarcodedIgReC!")
    except KeyboardInterrupt:
        log.info("\nBarcodedIgReC was interrupted!")
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

    log.info("Log was written to " + params.log_filename)


if __name__ == '__main__':
    main()
