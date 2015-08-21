#!/usr/bin/env python

import os
import shutil
import sys
import logging

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'
spades_src = os.path.join(home_directory, "src/python_pipeline/")

sys.path.append(spades_src)
import process_cfg
import support

#######################################################################################
#           Error messages
#######################################################################################

def ErrorMessagePrepareCfg(log):
    log.info("Probably you forgot to prepare IgRepertoireConstructor. Please follow instructions:")
    log.info("  (1) remove build/ directory")
    log.info("  (2) type command \'./prepare cfg\' to check all dependencies")
    log.info("  (3) type \'make\' to compile IgRepertoireConstrictor")
    log.info("  (4) rerun IgRepertoireConstrictor")

#######################################################################################
#           Binary routines
#######################################################################################

class IgBinaryConfig:
    def __init__(self):
        self.path_to_vj_aligner = 'build/release/bin/ig_kplus_vj_finder'
        self.run_vj_aligner = 'build/release/bin/./ig_kplus_vj_finder'
        self.path_to_trie_compressor = 'build/release/bin/ig_trie_compressor'
        self.run_trie_compressor = 'build/release/bin/./ig_trie_compressor'
        self.path_to_graph_constructor = 'build/release/bin/ig_swgraph_construct'
        self.run_graph_constructor = 'build/release/bin/./ig_swgraph_construct'
        self.path_to_consensus_constructor = 'build/release/bin/ig_final_alignment'
        self.run_consensus_constructor = 'build/release/bin/./ig_final_alignment'
        self.path_to_dsf = 'build/release/bin/dense_sgraph_finder'

    def CheckBinaries(self, log):
        if not os.path.exists(self.path_to_vj_aligner):
            log.info("ERROR: Binary file of VJFinder was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_trie_compressor):
            log.info("ERROR: Binary file of TrieCompressor was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_graph_constructor):
            log.info("ERROR: Binary file of GraphConstructor was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_dsf):
            log.info("ERROR: Binary file of DSF was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_consensus_constructor):
            log.info("ERROR: Binary file of ConsensusConstructor was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)

#######################################################################################
#           Phases
#######################################################################################
class PhaseNames:
    def __init__(self):
        self.__short_names = ['vj_alignment', 'trie_compressor', 'graph_constructor', 'dsf', 'consensus_constructor']
        self.__long_names = ['VJ Alignment', 'Trie Compressor', 'Graph Constructor', 'Dense Subgraph Finder', 'Consensus Constructor']

    def __iter__(self):
        for sname in self.__short_names:
            yield sname

    def __len__(self):
        return len(self.__short_names)

    def GetPhaseNameBy(self, index):
        return self.__short_names[index]

    def GetPhaseIndex(self, phase_name):
        for i in range(len(self)):
            if self.GetPhaseNameBy(i) == phase_name:
                return i
        return -1

    def PhaseIsVJAlignment(self, phase_name):
        return phase_name == 'vj_alignment'

    def PhaseIsTrieCompressor(self, phase_name):
        return phase_name == 'trie_compressor'

    def PhaseIsGraphConstructor(self, phase_name):
        return phase_name == 'graph_constructor'

    def PhaseIsDSF(self, phase_name):
        return phase_name == 'dsf'

    def PhaseIsConsensusConstructor(self, phase_name):
        return phase_name == 'consensus_constructor'

###########
class VJAlignmentPhase:
    def __init__(self, params, log):
        self.__params = params
        self.__log = log

    def PrintStartMessage(self):
        self.__log.info("==== VJ aligner starts")

    def Run(self):
        self.__params.vj_finder_output_dir = os.path.join(self.__params.output, "vj_finder")
        command_line = IgBinaryConfig().run_vj_aligner + " -i " + self.__params.reads + " -o " + \
                       self.__params.vj_finder_output_dir
        self.__log.info("VJ finder command line: " + command_line + "\n")
        support.sys_call(command_line, self.__log)

    def PrintOutputFiles(self):
        self.__log.info("Some output files!")

    def PrintFinishMessage(self):
        self.__log.info("==== VJ aligner finished\n")

###########
class TrieCompressionPhase:
    def __init__(self, params, log):
        self.__params = params
        self.__log = log

    def PrintStartMessage(self):
        self.__log.info("==== Trie Compressor starts")

    def Run(self):
        self.__log.info("Rrrrruning!")

    def PrintOutputFiles(self):
        self.__log.info("Some output files!")

    def PrintFinishMessage(self):
        self.__log.info("==== Trie Compressor finished\n")

###########
class GraphConstructionPhase:
    def __init__(self, params, log):
        self.__params = params
        self.__log = log

    def PrintStartMessage(self):
        self.__log.info("Hello, I am Graph Constructor!")

    def Run(self):
        self.__log.info("Rrrrruning!")

    def PrintOutputFiles(self):
        self.__log.info("Some output files!")

###########
class DSFPhase:
    def __init__(self, params, log):
        self.__params = params
        self.__log = log

    def PrintStartMessage(self):
        self.__log.info("Hello, I am Dense Subgraph Finder!")

    def Run(self):
        self.__log.info("Rrrrruning!")

    def PrintOutputFiles(self):
        self.__log.info("Some output files!")

###########
class ConsensusConstructionPhase:
    def __init__(self, params, log):
        self.__params = params
        self.__log = log

    def PrintStartMessage(self):
        self.__log.info("Hello, I am Consensus Constructor!")

    def Run(self):
        self.__log.info("Rrrrruning!")

    def PrintOutputFiles(self):
        self.__log.info("Some output files!")

###########
class PhaseFactory:
    def __init__(self, params, log):
        self.__entry_point = params.entry_point
        self.__params = params
        self.__log = log
        self.__phase_names = PhaseNames()

    def __CreatePhaseByName(self, phase_name):
        if self.__phase_names.PhaseIsVJAlignment(phase_name):
            return VJAlignmentPhase(self.__params, self.__log)
        elif self.__phase_names.PhaseIsTrieCompressor(phase_name):
            return TrieCompressionPhase(self.__params, self.__log)
        elif self.__phase_names.PhaseIsGraphConstructor(phase_name):
            return GraphConstructionPhase(self.__params, self.__log)
        elif self.__phase_names.PhaseIsDSF(phase_name):
            return DSFPhase(self.__params, self.__log)
        elif self.__phase_names.PhaseIsConsensusConstructor(phase_name):
            return ConsensusConstructionPhase(self.__params, self.__log)

    def CreatePhases(self):
        phase_list = list()
        first_phase_index = self.__phase_names .GetPhaseIndex(self.__entry_point)
        if first_phase_index == -1:
            self.__log.info("Incorrect name of entry-point")
            sys.exit(1)
        for i in range(first_phase_index, len(self.__phase_names )):
            phase_list.append(self.__CreatePhaseByName(self.__phase_names.GetPhaseNameBy(i)))
        return phase_list

############
class PhaseManager:
    def __init__(self, params, log):
        self.__params = params
        self.__log = log
        self.__phase_factory = PhaseFactory(params, log)
        self.__phases = self.__phase_factory.CreatePhases()

    def Run(self):
        for phase in self.__phases:
            phase.PrintStartMessage()
            phase.Run()
            phase.PrintOutputFiles()
            self.__log.info("\n============================================\n")

#######################################################################################
#           IO routines
#######################################################################################

def CreateLogger():
    # prepare log
    log = logging.getLogger('ig_repertoire_constructor')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    return log

def ParseCommandLineParams():
    # Parse commandline args
    from src.python_add.argparse_ext import ArgumentHiddenParser
    parser = ArgumentHiddenParser(description="IgRepertoireConstructor: an algorithm for construction of "
                                              "antibody repertoire from immunosequencing data",
                                  epilog="""
    In case you have troubles running IgRepertoireConstructor, you can write to igtools_support@googlegroups.com.
    Please provide us with ig_repertoire_constructor.log file from the output directory.
                                  """,
                                  add_help=False)

    req_args = parser.add_argument_group("Input")
    input_args = req_args.add_mutually_exclusive_group(required=True)
    input_args.add_argument("-s", "--reads",
                            type=str,
                            default="", # FIXME This is only for ace's version of python. Locally it works great w/o it
                            help="immunosequencing reads in FASTQ format")
    input_args.add_argument("--test",
                            action="store_const",
                            const="test_dataset/merged_reads.fastq",
                            dest="reads",
                            help="running of test dataset")

    out_args = parser.add_argument_group("Output")
    out_args.add_argument("-o", "--output",
                          type=str,
                          default="igrepcon_test",
                          help="output directory [default: \"%(default)s\"]")

    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument("-t", "--threads",
                               type=int,
                               default=16,
                               dest="num_threads",
                               help="thread number [default: %(default)d]")
    optional_args.add_argument("--tau",
                               type=int,
                               default=3,
                               dest="max_mismatches",
                               help="maximum allowed mismatches between identical error-prone reads [default: %(default)d]")
    optional_args.add_argument("-h", "--help",
                               action="help",
                               help="show this help message and exit")

    dev_args = parser.add_argument_group("_Developer arguments")
    dev_args.add_argument("-f", "--min-fillin",
                          type=float,
                          default=0.6,
                          help="_minimum edge fill-in of dense subgraphs [default: %(default)2.1f]")
    dev_args.add_argument('--entry-point',
                          type=str,
                          default=PhaseNames().GetPhaseNameBy(0),
                          help="_continue from the given stage [default: %(default)s]")
    ods_args = dev_args.add_mutually_exclusive_group(required=False)
    ods_args.add_argument("--help-hidden", "-H",
                          action="help_hidden",
                          help="_Show hidden help")
    parser.set_defaults(config_dir="configs",
                        config_file="config.info")
    return parser, parser.parse_args()

def CheckParamsCorrectness(parser, params, log):
    if not "output" in params or params.output == "":
        log.info("ERROR: Output directory (-o) was not specified\n")
        parser.print_help()
        sys.exit(-1)
    if not "reads" in params or params.reads == "":
        log.info("ERROR: Reads (-s) were not specified\n")
        parser.print_help()
        sys.exit(-1)
    if not os.path.exists(params.reads):
        log.info("ERROR: File with reads " + params.reads + " were not found\n")
        parser.print_help()
        sys.exit(-1)
    if not os.path.isabs(params.reads):
        params.reads = os.path.abspath(params.reads)

def PrepareOutputDir(params):
    if params.entry_point == "vj_alignment" and os.path.exists(params.output):
        shutil.rmtree(params.output)
    if not os.path.isdir(params.output):
        os.makedirs(params.output)

def PrintParams(params, log):
    log.info("IgRepertoireConstructor parameters:")
    log.info("  Input reads:\t\t\t" + params.reads)
    log.info("  Output directory:\t\t" + params.output)
    log.info("  Number of threads:\t\t" + str(params.num_threads))
    log.info("  Maximal number of mismatches:\t" + str(params.max_mismatches))
    log.info("  Entry point:\t\t\t" + params.entry_point)

def CreateFileLogger(params, log):
    params.log_filename = os.path.join(params.output, "ig_repertoire_constructor.log")
    if os.path.exists(params.log_filename):
        os.remove(params.log_filename)
    log_handler = logging.FileHandler(params.log_filename, mode='a')
    log.addHandler(log_handler)
    log.info("Log will be written to " + params.log_filename + "\n")

def PrintCommandLine(log):
    command_line = "Command line: " + " ".join(sys.argv)
    log.info("\n" + command_line + "\n")

#######################################################################################
#           Main
#######################################################################################
def main():
    binary_config = IgBinaryConfig()
    log = CreateLogger()
    binary_config.CheckBinaries(log)
    parser, params = ParseCommandLineParams()
    CheckParamsCorrectness(parser, params, log)
    PrepareOutputDir(params)
    CreateFileLogger(params, log)
    PrintParams(params, log)
    PrintCommandLine(log)

    # todo: add try-catch
    ig_repertoire_constructor = PhaseManager(params, log)
    ig_repertoire_constructor.Run()

if __name__ == '__main__':
    main()
