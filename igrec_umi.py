import os
import sys
from argparse import ArgumentParser

import igrec

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
spades_src = os.path.join(home_directory, "src/python_pipeline/")

sys.path.append(spades_src)
import support


def HelpAndReturn(log, parser, exit_code = 0):
    log.write = log.info
    parser.print_help(log)
    sys.exit(exit_code)


def ParseCommandLineParams(log):
    from argparse import ArgumentParser
    parser = ArgumentParser(description="IgReC: an algorithm for construction of antibody repertoire from immunosequencing data",
                            epilog="""
    In case you have troubles running IgReC, you can write to igtools_support@googlegroups.com.
    Please provide us with ig_repertoire_constructor.log file from the output directory.
                            """,
                            add_help=False)

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

    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument("-t", "--threads",
                               type=int,
                               default=16,
                               dest="num_threads",
                               help="Thread number [default: %(default)d]")
    optional_args.add_argument("--umi_graph_tau",
                               type=int,
                               default=1,
                               dest="umi_graph_tau",
                               help="Maximum allowed mismatches between UMIs, used for UMI correction"
                                    "[default: %(default)d]")
    optional_args.add_argument("--igrec_tau",
                               type=int,
                               default=2,
                               dest="igrec_tau",
                               help="Maximum allowed mismatches between UMI clusters"
                                    "[default: %(default)d]")
    optional_args.add_argument("-n", "--min-super-read-size",
                               type=int,
                               default=5,
                               dest="min_super_read_size",
                               help="Minimum super read size [default: %(default)d]")
    optional_args.set_defaults(compile = True)
    optional_args.add_argument("-p", "--no-compilation",
                               dest="no_compilation",
                               action="store_true",
                               help="Exclude c++ code compilation from the pipeline")
    optional_args.add_argument("-c", "--ignore_code",
                               dest="ignore_code_changes",
                               action="store_true",
                               help="Ignore code changes when checking stages depensences")

    vj_align_args = parser.add_argument_group("Algorithm arguments")
    vj_align_args.add_argument("-l", "--loci",
                               type=str,
                               dest="loci",
                               default="",
                               help="Loci: IGH, IGK, IGL, IG (all BCRs), TRA, TRB, TRG, TRD, TR (all TCRs) or all. Required")

    vj_align_args.add_argument("--organism",
                               type=str,
                               default="human",
                               dest="organism",
                               help="Organism (human and mouse only are supported for this moment) [default: %(default)s]")

    params = parser.parse_args()

    # process help
    if len(sys.argv) == 1:
        HelpAndReturn(log, parser)

    # Process pair reads
    if params.left_reads or params.right_reads:
        if not params.left_reads or not params.right_reads:
            log.info("ERROR: Both left (-1) and right (-2) paired-end reads should be specified\n")
            sys.exit(-1)
        params.single_reads = "%s/merged_reads.fastq" % params.output

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
                line = line.replace("%MIN_SUPER_NODE_SIZE", str(params.min_super_read_size))
                if '%' in line:
                    log.error("Not all template variables substituted in the makefile, update igrec_umi.py script, line #%d: '%s'" % (idx, line))
                    exit(1)
                file.write(line)

    @classmethod
    def Prepare(self, params, stage_template, log, stage_dest = None, makefile_name = MAKEFILE):
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
    _StagePrepare.Prepare(params, ".", log, makefile_name ="Makefile_vars")
    _StagePrepare.Prepare(params, "no_compilation" if params.no_compilation else "compilation", log, stage_dest ="compilation")
    if params.single_reads:
        _StagePrepare.Prepare(params, "vj_finder_input", log, stage_dest ="vj_finder")
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


def main():
    log = igrec.CreateLogger()
    parser, params = ParseCommandLineParams(log)
    CheckParamsCorrectness(parser, params, log)
    if not os.path.exists(params.output):
        os.makedirs(params.output)
    igrec.CreateFileLogger(params, log)
    igrec.PrintCommandLine(log)
    final_dir = InitMakeFiles(params, log)
    # We need freshly compiled version to get actual build info
    if not params.no_compilation:
        support.sys_call("make -C " + os.path.join(os.path.dirname(final_dir), "compilation"), log)
    from src.build_info.build_info import BuildInfo
    print "===================Build info==================="
    BuildInfo().Log(log)
    print "================================================"
    support.sys_call("make -C " + final_dir, log)


if __name__ == '__main__':
    main()
