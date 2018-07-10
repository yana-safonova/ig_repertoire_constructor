#!/usr/bin/env python2

import os
import shutil
import sys
import logging
from abc import ABCMeta, abstractmethod

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'

import dense_subgraph_finder

import support

#######################################################################################
#           Error messages
#######################################################################################

def ErrorMessagePrepareCfg(log):
    log.info("Probably you forgot to prepare IgReC. Please follow instructions:")
    log.info("  (1) remove build/ directory")
    log.info("  (2) type command \'./prepare cfg\' to check all dependencies")
    log.info("  (3) type \'make\' to compile IgReC")
    log.info("  (4) rerun IgReC")

def SupportInfo(log):
    log.info("\nIn case you have troubles running IgReC, "
             "you can write to igtools_support@googlegroups.com.")
    log.info("Please provide us with igrec.log file from the output directory.")

#######################################################################################
#           Binary routines
#######################################################################################
class IgRepConConfig:
    def __initBinaryPaths(self):
        self.path_to_pair_reads_merger = os.path.join(home_directory, 'build/release/bin/paired_read_merger')
        self.run_pair_reads_merger = os.path.join(home_directory, 'build/release/bin/./paired_read_merger')
        self.path_to_vj_aligner = os.path.join(home_directory, 'build/release/bin/vj_finder') #ig_kplus_vj_finder')
        self.run_vj_aligner = os.path.join(home_directory, 'build/release/bin/./vj_finder') #ig_kplus_vj_finder')
        self.path_to_trie_compressor = os.path.join(home_directory, 'build/release/bin/ig_trie_compressor')
        self.run_trie_compressor = os.path.join(home_directory, 'build/release/bin/./ig_trie_compressor')
        self.run_fake_trie_compressor = os.path.join(home_directory, 'py/ig_fake_trie_compressor.py')
        self.path_to_graph_constructor = os.path.join(home_directory, 'build/release/bin/ig_swgraph_construct')
        self.run_graph_constructor = os.path.join(home_directory, 'build/release/bin/./ig_swgraph_construct')
        self.path_to_consensus_constructor = os.path.join(home_directory, 'build/release/bin/ig_component_splitter')
        self.run_consensus_constructor = os.path.join(home_directory, 'build/release/bin/./ig_component_splitter')
        self.run_rcm_recoverer = os.path.join(home_directory, 'py/rcm_recoverer.py')
        self.run_compress_equal_clusters = os.path.join(home_directory, 'py/ig_compress_equal_clusters.py')
        self.run_report_supernodes = os.path.join(home_directory, 'py/ig_report_supernodes.py')
        self.run_triecmp_to_repertoire = os.path.join(home_directory, 'py/ig_triecmp_to_repertoire.py')
        self.run_divan = os.path.join(home_directory, 'diversity_analyzer.py')
        self.run_to_vidjil = os.path.join(home_directory, 'py/to_vidjil.py')
        self.path_to_dsf = os.path.join(home_directory, 'build/release/bin/dense_sgraph_finder')
        self.path_to_divan = os.path.join(home_directory, 'build/release/bin/cdr_labeler')
        self.path_to_germline = os.path.join(home_directory, "data/germline")

    def __init__(self):
        self.__initBinaryPaths()

    def CheckBinaries(self, log):
        if not os.path.exists(self.path_to_pair_reads_merger):
            log.info("ERROR: Binary file of " + PairReadMergerPhase.GetLongName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_vj_aligner):
            log.info("ERROR: Binary file of " + VJAlignmentPhase.GetLongName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_trie_compressor):
            log.info("ERROR: Binary file of " + TrieCompressionPhase.GetLongName() + " (" + self.path_to_trie_compressor +") was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.run_report_supernodes):
            log.info("ERROR: Binary file of " + TrieCompressionPhase.GetLongName() +  " (" + self.run_report_supernodes + ") was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.run_triecmp_to_repertoire):
            log.info("ERROR: Binary file of " + TrieCompressionPhase.GetLongName() + " (" + self.run_triecmp_to_repertoire + ") was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_graph_constructor):
            log.info("ERROR: Binary file of " + GraphConstructionPhase.GetLongName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_dsf):
            log.info("ERROR: Binary file of " + DSFPhase.GetLongName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_consensus_constructor):
            log.info("ERROR: Binary file of " + ConsensusConstructionPhase.GetLongName() + " (" + self.path_to_consensus_constructor + ") was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.run_rcm_recoverer):
            log.info("ERROR: Binary file of " + ConsensusConstructionPhase.GetLongName() + " (" + self.run_rcm_recoverer + ") was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.run_compress_equal_clusters):
            log.info("ERROR: Binary file of " + CompressEqualClustersPhase.GetLongName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.run_divan) or not os.path.exists(self.path_to_divan):
            log.info("ERROR: Binary file of " + DiversityAnalyzerPhase.GetLongName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.run_to_vidjil):
            log.info("ERROR: Binary file of " + DiversityAnalyzerPhase.GetLongName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)


class IgRepConIO:
    def __initVJFinderOutput(self, output_dir):
        self.vj_finder_output = os.path.join(output_dir, "vj_finder")
        self.cropped_reads = os.path.join(self.vj_finder_output, "cleaned_reads.fa")
        self.bad_reads = os.path.join(self.vj_finder_output, "filtered_reads.fa")
        self.vj_alignment_info = os.path.join(self.vj_finder_output, "alignment_info.csv")

    def __initCompressorOutput(self, output_dir):
        self.compressed_reads = os.path.join(output_dir, "compressed_reads.fa")
        self.map_file = os.path.join(output_dir, "cleaned_compressed_map.txt")
        self.supernodes_file = os.path.join(output_dir, "super_reads.fa")
        self.supernode_repertoire = os.path.join(output_dir, "supernode_repertoire.fa")
        self.supernode_rcm = os.path.join(output_dir, "supernode_repertoire.rcm")

    def __initDSFOutput(self, output_dir):
        self.dsf_output = os.path.join(output_dir, "dense_sgraph_finder")
        self.dense_sgraph_decomposition = os.path.join(self.dsf_output, 'dense_subgraphs.txt')

    def __initFinalOutput(self, output_dir):
        self.uncompressed_final_clusters_fa = os.path.join(output_dir, 'final_repertoire_uncompressed.fa')
        self.uncompressed_final_rcm = os.path.join(output_dir, 'final_repertoire_uncompressed.rcm')

    def __initCompressEqualClusters(self, output_dir):
        self.tmp_compressed_clusters_fa = os.path.join(output_dir, 'tmp_compressed_clusters.fa')
        self.tmp_compressed_clusters_map = os.path.join(output_dir, 'tmp_compressed_clusters.map')
        self.compressed_final_clusters_fa = os.path.join(output_dir, 'final_repertoire.fa')
        self.compressed_final_rcm = os.path.join(output_dir, 'final_repertoire.rcm')

    def __initDiversityAnalyzer(self, output_dir):
        self.divan_output = os.path.join(output_dir, 'divan')
        self.vidjil_output = os.path.join(output_dir, 'igrec.vidjil')
        self.output = output_dir
        self.divan_feature_file = os.path.join(self.divan_output, 'cdr_details.txt')

    def __init__(self, output_dir, log):
        self.__log = log
        self.__initVJFinderOutput(output_dir)
        self.__initCompressorOutput(output_dir)
        self.sw_graph = os.path.join(output_dir, "sw.graph")
        self.__initDSFOutput(output_dir)
        self.__initFinalOutput(output_dir)
        self.final_stripped_clusters_fa = os.path.join(output_dir, 'final_repertoire_large.fa')
        self.__initCompressEqualClusters(output_dir)
        self.__initDiversityAnalyzer(output_dir)


#######################################################################################
#           Phases
#######################################################################################
class Phase:
    __metaclass__ = ABCMeta

    def __init__(self, log):
        self._log = log

    @classmethod
    def GetLongName(cls):
        raise NotImplementedError()

    def Run(self):
        self._PrintStartMessage()
        self._CheckFilesExistence(self._GetInputFiles())
        self._Run()
        self._CheckFilesExistence(self._GetOutputFiles())
        self._PrintOutputFiles()
        self._PrintFinishMessage()

    def _PrintStartMessage(self):
        self._log.info("==== " + self.GetLongName() + " starts\n")

    def _PrintFinishMessage(self):
        self._log.info("\n==== " + self.GetLongName() + " finished")

    @abstractmethod
    def _GetInputFiles(self):
        pass

    @abstractmethod
    def _GetOutputFiles(self):
        pass

    def _CheckFilesExistence(self, files_with_descr):
        for (file_path, description) in files_with_descr:
            if not os.path.exists(file_path):
                self._log.error('Could not find the file with ' + description + ' at ' + file_path)
                SupportInfo(self._log)
                sys.exit(1)

    def _PrintOutputFiles(self):
        self._log.info("\nOutput files: ")
        for (file, description) in self._GetOutputFiles():
            self._log.info('  * ' + description + ' can be found at ' + file)

    @abstractmethod
    def _Run(self):
        pass


###########
class PairReadMergerPhase(Phase):
    def __init__(self, params, log):
        super(PairReadMergerPhase, self).__init__(log)
        self.__params = params

    @classmethod
    def GetLongName(cls):
        return 'Pair reads merging'

    def _GetInputFiles(self):
        return (
            (self.__params.left_reads, 'input left reads'),
            (self.__params.right_reads, 'input right reads'),
        )

    def _GetOutputFiles(self):
        return (
            (self.__params.single_reads, 'input reads'),
        )

    def _Run(self):
        command_line = "%s %s %s %s" % (IgRepConConfig().run_pair_reads_merger,
                                        self.__params.left_reads,
                                        self.__params.right_reads,
                                        self.__params.single_reads)
        cpuprofile = self.__params.output + "/pair_read_merger_prof.out" if self.__params.profile else None
        support.sys_call_ex(command_line, self._log, cpuprofile=cpuprofile)


###########
class VJAlignmentPhase(Phase):
    def __init__(self, params, log):
        super(VJAlignmentPhase, self).__init__(log)
        self.__params = params

    @classmethod
    def GetLongName(cls):
        return 'VJ Alignment'

    def _GetInputFiles(self):
        return (
            (self.__params.single_reads, 'input reads'),
        )

    def _GetOutputFiles(self):
        output_files = [
            (self.__params.io.cropped_reads, 'cleaned Ig-Seq reads'),
        ]
        if not self.__params.no_alignment:
            output_files.extend([
                (self.__params.io.bad_reads, 'contaminated (not Ig-Seq) reads'),
                (self.__params.io.vj_alignment_info, 'VJ alignment output'),
            ])
        return output_files

    def _Run(self):
        if not self.__params.no_alignment:
            self.__params.vj_finder_output = os.path.join(self.__params.output, "vj_finder")
            command_line = os.path.abspath(IgRepConConfig().run_vj_aligner) + \
                " -i " + os.path.abspath(self.__params.single_reads) + \
                " -o " + os.path.abspath(self.__params.io.vj_finder_output) + \
                " --db-directory " + os.path.abspath(IgRepConConfig().path_to_germline) + \
                " -t " + str(self.__params.num_threads) + \
                " --loci " + self.__params.loci + \
                " --organism " + self.__params.organism
            if self.__params.no_pseudogenes:
                command_line += " --pseudogenes=off"
            else:
                command_line += " --pseudogenes=on"

            cpuprofile = self.__params.output + "/vjf_prof.out" if self.__params.profile else None

            cwd = os.getcwd()
            os.chdir(home_directory)
            support.sys_call_ex(command_line, self._log, cpuprofile=cpuprofile)
            os.chdir(cwd)
        else:
            self._log.info("VJ Finder stage skipped")
            self.__params.io.cropped_reads = self.__params.single_reads


###########
class TrieCompressionPhase(Phase):
    def __init__(self, params, log):
        super(TrieCompressionPhase, self).__init__(log)
        self.__params = params

    @classmethod
    def GetLongName(cls):
        return 'Trie Compressor'

    def _GetInputFiles(self):
        return (
            (self.__params.io.cropped_reads, 'cleaned Ig-Seq reads'),
        )

    def _GetOutputFiles(self):
        return (
            (self.__params.io.compressed_reads, 'compressed reads'),
            (self.__params.io.supernodes_file, 'super reads'),
            (self.__params.io.map_file, 'map from cleaned reads to compressed reads'),
        )

    def _Run(self):
        command_line = IgRepConConfig().run_trie_compressor + " -i " + self.__params.io.cropped_reads + \
                    " -o " + self.__params.io.compressed_reads + " -m " + self.__params.io.map_file + " -Toff"
        cpuprofile = self.__params.output + "/trie_compressor_prof.out" if self.__params.profile else None
        support.sys_call_ex(command_line, self._log, cpuprofile=cpuprofile)

        command_line = IgRepConConfig().run_triecmp_to_repertoire + " -i " + self.__params.io.cropped_reads + \
                       " -c " + self.__params.io.compressed_reads + " -m " + self.__params.io.map_file + \
                       " -r " + self.__params.io.supernode_repertoire + " -R " + self.__params.io.supernode_rcm
        support.sys_call(command_line, self._log)
        command_line = "%s %s %s --limit=%d" % (IgRepConConfig().run_report_supernodes,
                                                self.__params.io.supernode_repertoire,
                                                self.__params.io.supernodes_file,
                                                self.__params.min_cluster_size)
        support.sys_call(command_line, self._log)

        if not self.__params.equal_compression:
            command_line = IgRepConConfig().run_fake_trie_compressor + " -i " + self.__params.io.cropped_reads + \
                        " -o " + self.__params.io.compressed_reads + " -m " + self.__params.io.map_file
            support.sys_call(command_line, self._log)


###########
class GraphConstructionPhase(Phase):
    def __init__(self, params, log):
        super(GraphConstructionPhase, self).__init__(log)
        self.__params = params

    @classmethod
    def GetLongName(cls):
        return 'Graph Constructor'

    def _GetInputFiles(self):
        return (
            (self.__params.io.compressed_reads, 'compressed reads'),
        )

    def _GetOutputFiles(self):
        return (
            (self.__params.io.sw_graph, 'Smith-Waterman graph'),
        )

    def _Run(self):
        command_line = IgRepConConfig().run_graph_constructor + " -i " + self.__params.io.compressed_reads + \
                       " -o " + self.__params.io.sw_graph + " -t " + str(self.__params.num_threads) + \
                       " --tau=" + str(self.__params.max_mismatches) + " -A" + " -Toff"
        cpuprofile = self.__params.output + "/graph_constructor_prof.out" if self.__params.profile else None
        support.sys_call_ex(command_line, self._log, cpuprofile=cpuprofile)


###########
class DSFPhase(Phase):
    def __init__(self, params, log):
        super(DSFPhase, self).__init__(log)
        self.__params = params

    @classmethod
    def GetLongName(cls):
        return 'Dense Subgraph Finder'

    def _GetInputFiles(self):
        return (
            (self.__params.io.sw_graph, 'Smith-Waterman graph'),
        )

    def _GetOutputFiles(self):
        return (
            (self.__params.io.dense_sgraph_decomposition, 'dense subgraph decomposition'),
        )

    def __GetDSFParams(self):
        dsf_params = ['-g', self.__params.io.sw_graph,
                      '-o', self.__params.io.dsf_output,
                      '-t', str(self.__params.num_threads),
                      '-n', str(self.__params.min_snode_size),
                      '-f', str(self.__params.min_fillin)]
        if self.__params.create_trivial_decomposition:
            dsf_params.append('--create-triv-dec')
        if self.__params.save_aux_files:
            dsf_params.append('--save-aux-files')
        return dsf_params

    def _Run(self):
        dense_subgraph_finder.main(self.__GetDSFParams(), self.__params.log_filename)


###########
class ConsensusConstructionPhase(Phase):
    def __init__(self, params, log):
        super(ConsensusConstructionPhase, self).__init__(log)
        self.__params = params

    @classmethod
    def GetLongName(cls):
        return 'Consensus Constructor'

    def _GetInputFiles(self):
        return (
            (self.__params.io.compressed_reads, 'compressed reads'),
            (self.__params.io.dense_sgraph_decomposition, 'dense subgraph decomposition'),
            (self.__params.io.cropped_reads, 'cleaned Ig-Seq'),
            (self.__params.io.map_file, 'map from cleaned reads to compressed reads'),
        )

    def _GetOutputFiles(self):
        return (
            (self.__params.io.uncompressed_final_clusters_fa, 'antibody clusters of uncompressed final repertoire'),
            (self.__params.io.uncompressed_final_rcm, 'read-cluster map of uncompressed final repertoire'),
        )

    def _Run(self):
        command_line = "%s -i %s -c %s -q %s -o %s" % (IgRepConConfig().run_rcm_recoverer,
                                                       self.__params.io.cropped_reads,
                                                       self.__params.io.map_file,
                                                       self.__params.io.dense_sgraph_decomposition,
                                                       self.__params.io.uncompressed_final_rcm)
        support.sys_call(command_line, self._log)
        command_line = IgRepConConfig().run_consensus_constructor + \
                       " -i " + self.__params.io.cropped_reads + \
                       " -R " + self.__params.io.uncompressed_final_rcm + \
                       " -M " + self.__params.io.uncompressed_final_rcm + \
                       " -o " + self.__params.io.uncompressed_final_clusters_fa + \
                       " -t " + str(self.__params.num_threads) + \
                       " -D " + str(self.__params.discard) + \
                       " --max-votes " + str(self.__params.max_votes)
        cpuprofile = self.__params.output + "/consensus_constructor_prof.out" if self.__params.profile else None
        support.sys_call_ex(command_line, self._log, cpuprofile=cpuprofile)


class CompressEqualClustersPhase(Phase):
    def __init__(self, params, log):
        super(CompressEqualClustersPhase, self).__init__(log)
        self.__params = params

    @classmethod
    def GetLongName(cls):
        return 'Compress Equal Final Clusters'

    def _GetInputFiles(self):
        return (
            (self.__params.io.uncompressed_final_clusters_fa, 'antibody clusters of uncompressed final repertoire'),
            (self.__params.io.uncompressed_final_rcm, 'read-cluster map of uncompressed final repertoire'),
        )

    def _GetOutputFiles(self):
        return (
            (self.__params.io.compressed_final_clusters_fa, 'compressed antibody clusters of final repertoire'),
            (self.__params.io.compressed_final_rcm, 'read-cluster map of compressed final repertoire'),
        )

    def _Run(self):
        command_line = "%s %s %s -T %s -m %s -r %s -R %s" % (IgRepConConfig().run_compress_equal_clusters,
                                                             self.__params.io.uncompressed_final_clusters_fa,
                                                             self.__params.io.compressed_final_clusters_fa,
                                                             self.__params.io.tmp_compressed_clusters_fa,
                                                             self.__params.io.tmp_compressed_clusters_map,
                                                             self.__params.io.uncompressed_final_rcm,
                                                             self.__params.io.compressed_final_rcm)
        support.sys_call(command_line, self._log)


class RemoveLowAbundanceReadsPhase(Phase):
    def __init__(self, params, log):
        super(RemoveLowAbundanceReadsPhase, self).__init__(log)
        self.__params = params

    @classmethod
    def GetLongName(cls):
        return 'Low Abundant Clusters Remover'

    def _GetInputFiles(self):
        return (
            (self.__params.io.compressed_final_clusters_fa, 'compressed antibody clusters of final repertoire'),
        )

    def _GetOutputFiles(self):
        return (
            (self.__params.io.final_stripped_clusters_fa, 'highly abundant antibody clusters of final repertoire'),
        )

    def _Run(self):
        command_line = "%s %s %s --limit=%d" % (IgRepConConfig().run_report_supernodes,
                                                self.__params.io.compressed_final_clusters_fa,
                                                self.__params.io.final_stripped_clusters_fa,
                                                self.__params.min_cluster_size)
        support.sys_call(command_line, self._log)


class DiversityAnalyzerPhase(Phase):
    def __init__(self, params, log):
        super(DiversityAnalyzerPhase, self).__init__(log)
        self.__params = params

    @classmethod
    def GetLongName(cls):
        return 'IgDiversityAnalyzer'

    def _GetInputFiles(self):
        return (
            (self.__params.io.compressed_final_clusters_fa, 'compressed antibody clusters of final repertoire'),
        )

    def _GetOutputFiles(self):
        return (
            (self.__params.io.divan_feature_file, 'repertoire sequence features'),
        )

    def _Run(self):
        command_line = "%s -i %s -t %d -o %s -l %s --org %s" % (
            IgRepConConfig().run_divan,
            self.__params.io.compressed_final_clusters_fa,
            self.__params.num_threads,
            self.__params.io.divan_output,
            self.__params.loci,
            self.__params.organism
        )
        support.sys_call(command_line, self._log)
        command_line = "%s -i %s -o %s --initial-reads %s" % (
            IgRepConConfig().run_to_vidjil,
            self.__params.io.output,
            self.__params.io.vidjil_output,
            self.__params.single_reads)
        support.sys_call(command_line, self._log)


###########
class PhaseFactory:

    __phase_order = (
        (PairReadMergerPhase, 'pair_reads_merger'),
        (VJAlignmentPhase, 'vj_alignment'),
        (TrieCompressionPhase, 'trie_compressor'),
        (GraphConstructionPhase, 'graph_constructor'),
        (DSFPhase, 'dsf'),
        (ConsensusConstructionPhase, 'consensus_constructor'),
        (CompressEqualClustersPhase, 'compress_equal_clusters'),
        (RemoveLowAbundanceReadsPhase, 'remove_low_abundance_reads'),
        (DiversityAnalyzerPhase, 'diversity_analyzer'),
    )

    def __init__(self, params, log):
        self.__entry_point = params.entry_point
        if self.__entry_point is None:
            if params.left_reads:
                self.__entry_point = self.__phase_order[0][1]
            else:
                self.__entry_point = self.__phase_order[1][1]
        self.__params = params
        self.__log = log

    def __CreatePhaseByName(self, phase_name):
        return next(phase for phase, phase_id in self.__phase_order if phase_id == phase_name)()

    def CreatePhases(self):
        phase_ids = [phase_id for phase, phase_id in self.__phase_order]
        if self.__entry_point in phase_ids:
            first_phase_index = next(idx for idx, phase_id in enumerate(phase_ids) if phase_id == self.__entry_point)
        else:
            self.__log.info("Incorrect name of entry-point")
            sys.exit(1)
        return [phase(self.__params, self.__log) for phase, _ in self.__phase_order[first_phase_index :]]


############
class PhaseManager:
    def __init__(self, phase_factory, params, log):
        self.__log = log
        self.__phases = phase_factory.CreatePhases()

    def __PrintPhaseDelimeter(self):
        self.__log.info("\n============================================\n")

    def Run(self):
        first = True
        for phase in self.__phases:
            if not first:
                self.__PrintPhaseDelimeter()
            first = False
            phase.Run()


#######################################################################################
#           IO routines
#######################################################################################

def CreateLogger():
    log = logging.getLogger('igrec')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.INFO)
    log.addHandler(console)
    return log

def HelpString():
    return "Usage: igrec.py (-s FILENAME | -1 FILENAME -2 FILENAME | --test)\n" +\
    "                (-o OUTPUT_DIR) (-l LOCI)\n" +\
    "                [-t / --threads INT]\n" +\
    "                [--organism ORGANISM] [--no-pseudogenes]\n" +\
    "                [--tau INT] [--min-sread-size INT] [--min-cluster-size INT]\n" +\
    "                [-h]\n\n" +\
    "IgReC: an algorithm for construction of antibody repertoire from immunosequencing data\n\n" +\
    "Input arguments:\n" +\
    "  -s\t\t\t\tFILENAME\t\tSingle reads in FASTQ format\n" +\
    "  -1\t\t\t\tFILENAME\t\tLeft paired-end reads in FASTQ format\n" +\
    "  -2\t\t\t\tFILENAME\t\tRight paired-end reads in FASTQ format\n" +\
    "  --test\t\t\t\t\t\tRunning of test dataset\n\n" +\
    "Output arguments:\n" +\
    "  -o / --output\t\t\tOUTPUT_DIR\t\tOutput directory. Required\n\n" +\
    "Running arguments:\n" +\
    "  -t / --threads\t\tINT\t\t\tThread number [default: 16]\n" +\
    "  -h / --help\t\t\t\t\t\tShowing help message and exit\n\n" +\
    "Alignment arguments:\n" +\
    "  --no-alignment\t\t\t\t\tDo not provide any alignment and filtering\n" +\
    "  -l / --loci\t\t\tLOCI\t\t\tLoci: IGH, IGK, IGL, IG (all BCRs), TRA, TRB, TRG, TRD, TR (all TCRs) or all. Required\n" +\
    "  --organism\t\t\tORGANISM\t\tOrganism: human, mouse, pig, rabbit, rat, rhesus_monkey are available [default: human]\n" +\
    "  --no-pseudogenes\t\t\t\t\tDisabling using pseudogenes along with normal gene segments for VJ alignment [default: False]\n\n" +\
    "Algorithm arguments:\n" +\
    "  --tau\t\t\t\tINT\t\t\tMaximum allowed mismatches between identical error-prone reads [default: 4]\n" +\
    "  --n / --min-sread-size\tINT\t\t\tMinimum size of super reads [default: 5]\n" +\
    "  --min-cluster-size\t\tINT\t\t\tMinimum size of antibody cluster using for output of large antibody clusters [default: 5]\n\n" +\
    "In case you have troubles running IgReC, you can write to igtools_support@googlegroups.com.\n" +\
    "Please provide us with igrec.log file from the output directory."


def ParseCommandLineParams(log):
    import argparse
    parser = argparse.ArgumentParser(description="IgReC: an algorithm for construction of "
                                     "antibody repertoire from immunosequencing data",
                                     epilog="""
    In case you have troubles running IgReC, you can write to igtools_support@googlegroups.com.
    Please provide us with ig_repertoire_constructor.log file from the output directory.
                                     """,
                                     add_help=False)

    class ActionTest(argparse.Action):
        def __init__(self, option_strings, dest, nargs=None, **kwargs):
            super(ActionTest, self).__init__(option_strings, dest, nargs=0, **kwargs)

        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, "single_reads", os.path.join(home_directory, "test_dataset/merged_reads.fastq"))
            setattr(namespace, "loci", "IG")
            setattr(namespace, "output", "igrec_test")

    req_args = parser.add_argument_group("Input")
    input_args = req_args.add_mutually_exclusive_group(required=False)
    input_args.add_argument("-s",
                            dest="single_reads",
                            type=str,
                            default="",  # FIXME This is only for ace's version of python. Locally it works great w/o it
                            help="Single reads in FASTQ format")

    input_args.add_argument("--test",
                            action=ActionTest,
                            default="",
                            help="Running of test dataset")

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
    optional_args.add_argument("--tau",
                               type=int,
                               default=4,
                               dest="max_mismatches",
                               help="Maximum allowed mismatches between identical error-prone reads "
                                    "[default: %(default)d]")
    optional_args.add_argument("--min-cluster-size",
                               type=int,
                               dest="min_cluster_size",
                               default=5,
                               help="Minimal size of antibody cluster using for output of large antibody clusters [default: %(default)d]")
    optional_args.add_argument("-n", "--min-sread-size",
                               type=int,
                               default=5,
                               dest="min_snode_size",
                               help="Minimum super read size [default: %(default)d]")

    optional_args.add_argument("--profile",
                               action="store_true",
                               default=False,
                               help="Enable CPU profiling")
    optional_args.add_argument("-h", "--help",
                               action="store_const",
                               const=True,
                               dest="show_help",
                               help="Showing help message and exit")

    vj_align_args = parser.add_argument_group("Algorithm arguments")
    vj_align_args.add_argument("-l", "--loci",
                               type=str,
                               dest="loci",
                               default="",
                               help="Loci: IGH, IGK, IGL, IG (all BCRs), TRA, TRB, TRG, TRD, TR (all TCRs) or all. Required")

    vj_align_args.add_argument("--no-pseudogenes",
                               action="store_const",
                               const=True,
                               dest="no_pseudogenes",
                               help="Do not use pseudogenes along with normal gene segments for VJ alignment [default: False]")

    vj_align_args.add_argument("--organism",
                               type=str,
                               default="human",
                               dest="organism",
                               help="Organism (human and mouse only are supported for this moment) [default: %(default)s]")

    vj_align_args.add_argument("--no-alignment",
                               action="store_true",
                               help="Do not provide any alignment and filtering")

    dev_args = parser.add_argument_group("Developer arguments")
    dev_args.add_argument("-f", "--min-fillin",
                          type=float,
                          default=0.6,
                          help="Minimum edge fill-in of dense subgraphs [default: %(default)2.1f]")
    dev_args.add_argument('--entry-point',
                          type=str,
                          default=None,
                          help="Continue from the given stage [default: %(default)s]")
    dev_args.add_argument("--create-triv-dec",
                          action="store_const",
                          const=True,
                          dest="create_trivial_decomposition",
                          help='Creating decomposition according to connected components [default: False]')
    dev_args.add_argument("--save-aux-files",
                          action="store_const",
                          const=True,
                          dest="save_aux_files",
                          help="Saving auxiliary files: subgraphs in GRAPH format and their decompositions "
                          "[default: False]")
    dev_args.add_argument("--debug",
                          action="store_const",
                          const=True,
                          dest="debug_mode",
                          help="Save auxiliary files [default: False]")
    dev_args.add_argument("-V", "--max-votes",
                          type=int,
                          default=10005000,
                          help="Maximun secondary votes threshold [default: %(default)d]")
    dev_args.add_argument("-D", "--discard",
                          action="store_true",
                          dest="discard",
                          help="Discard seconary vote clusters")
    dev_args.add_argument("--no-discard",
                          action="store_false",
                          help="Do not discard seconary vote clusters (default)")
    parser.set_defaults(discard=True)
    # TODO Add it into the help
    dev_args.add_argument("--no-equal-compression",
                          action="store_false",
                          dest="equal_compression",
                          help="Disable equal read compression before graph construction")
    dev_args.add_argument("--equal-compression",
                          action="store_true",
                          dest="equal_compression",
                          help="Enable equal read compression before graph construction (default)")
    parser.set_defaults(equal_compression=True)

    ods_args = dev_args.add_mutually_exclusive_group(required=False)
    ods_args.add_argument("--help-hidden", "-H",
                          action="help",
                          help="Show hidden help")
    parser.set_defaults(config_dir="configs",
                        config_file="config.info")
    params = parser.parse_args()

    # process help
    if params.show_help or len(sys.argv) == 1:
        log.info(HelpString())
        sys.exit(0)

    # Process pair reads
    if params.left_reads or params.right_reads:
        if not params.left_reads or not params.right_reads:
            log.info("ERROR: Both left (-1) and right (-2) paired-end reads should be specified\n")
            sys.exit(-1)
        params.single_reads = "%s/merged_reads.fastq" % params.output

    return parser, params


def EnsureAbsPath(s):
    if not os.path.isabs(s):
        s = os.path.abspath(s)
    return s


def CheckGeneralParamsCorrectness(parser, params, log):
    if not "output" in params or params.output == "":
        log.info("ERROR: Output directory (-o) was not specified\n")
        HelpString()
        sys.exit(-1)
    if not params.no_alignment and ("loci" not in params or params.loci == ""):
        log.info("ERROR: Immunological loci (-l) was not specified\n")
        HelpString()
        sys.exit(1)


def CheckSingleReadsCorrectness(parser, params, log):
    if not "single_reads" in params or params.single_reads == "":
        log.info("ERROR: Single reads (-s) were not specified\n")
        HelpString()
        sys.exit(-1)
    if not os.path.exists(params.single_reads):
        log.info("ERROR: File with single reads " + params.single_reads + " were not found\n")
        HelpString()
        sys.exit(-1)
    params.single_reads = EnsureAbsPath(params.single_reads)


def CheckPairedReadsCorrectness(parser, params, log):
    if not "left_reads" in params or params.left_reads == "":
        log.info("ERROR: Left reads (-1) were not specified\n")
        HelpString()
        sys.exit(-1)
    if not "right_reads" in params or params.right_reads == "":
        log.info("ERROR: Right reads (-2) were not specified\n")
        HelpString()
        sys.exit(-1)
    if not os.path.exists(params.left_reads):
        log.info("ERROR: File with left reads " + params.left_reads + " were not found\n")
        HelpString()
        sys.exit(-1)
    if not os.path.exists(params.right_reads):
        log.info("ERROR: File with right reads " + params.right_reads + " were not found\n")
        HelpString()
        sys.exit(-1)
    params.left_reads = EnsureAbsPath(params.left_reads)
    params.right_reads = EnsureAbsPath(params.right_reads)


def PrepareOutputDir(params):
    if params.entry_point == "vj_alignment" and os.path.exists(params.output):
        shutil.rmtree(params.output)
    if not os.path.isdir(params.output):
        os.makedirs(params.output)


def PrintParams(params, log):
    log.info("IgReC parameters:")
    log.info("  Input reads:\t\t\t" + params.single_reads)
    log.info("  Output directory:\t\t" + params.output)
    log.info("  Number of threads:\t\t" + str(params.num_threads))
    log.info("  Maximal number of mismatches:\t" + str(params.max_mismatches))
    log.info("  Entry point:\t\t\t" + params.entry_point if params.entry_point is not None else "start")


def CreateFileLogger(params, log):
    params.log_filename = os.path.join(params.output, "igrec.log")
    if os.path.exists(params.log_filename):
        os.remove(params.log_filename)
    log_handler = logging.FileHandler(params.log_filename, mode='a')
    log_handler.setLevel(logging.DEBUG)
    log.addHandler(log_handler)
    log.info("Log will be written to " + params.log_filename + "\n")


def PrintCommandLine(log):
    command_line = "Command line: " + " ".join(sys.argv)
    log.info("\n" + command_line + "\n")


def RemoveAuxFiles(params):
    if params.debug_mode:
        return
    if os.path.exists(params.io.map_file):
        os.remove(params.io.map_file)
    if os.path.exists(params.io.compressed_reads):
        os.remove(params.io.compressed_reads)
    if os.path.exists(params.io.sw_graph):
        os.remove(params.io.sw_graph)
    if os.path.exists(params.io.dsf_output) and not params.save_aux_files:
        shutil.rmtree(params.io.dsf_output)
    if os.path.exists(params.io.uncompressed_final_clusters_fa):
        os.remove(params.io.uncompressed_final_clusters_fa)
    if os.path.exists(params.io.uncompressed_final_rcm):
        os.remove(params.io.uncompressed_final_rcm)
    #if os.path.exists(params.io.merged_reads)


def PrintOutputFiles(params, log):
    log.info("\nIgReC output:")
    if os.path.exists(params.io.cropped_reads):
        log.info("  * Cleaned Ig-Seq reads were written to " + params.io.cropped_reads)
    if os.path.exists(params.io.bad_reads):
        log.info("  * Contaminated (not Ig-Seq) reads were written to " + params.io.bad_reads)
    if os.path.exists(params.io.vj_alignment_info):
        log.info("  * VJ alignment output was written to " + params.io.vj_alignment_info)
    if os.path.exists(params.io.supernodes_file):
        log.info("  * Super reads were written to " + params.io.supernodes_file)
    if os.path.exists(params.io.compressed_final_clusters_fa):
        log.info("  * Antibody clusters of final repertoire were written to " + params.io.compressed_final_clusters_fa)
    if os.path.exists(params.io.compressed_final_rcm):
        log.info("  * Read-cluster map of final repertoire was written to " + params.io.compressed_final_rcm)
    if os.path.exists(params.io.final_stripped_clusters_fa):
        log.info("  * Highly abundant antibody clusters of final repertoire were written to " + params.io.final_stripped_clusters_fa)


def LogInfo(log):
    import sys
    sys.path.append(home_directory + "/py")
    import build_info
    from datetime import datetime

    class LogInfoToDebugAdapter:
        def __init__(self, log):
            self.log__ = log
        def info(self, *args, **kwargs):
            self.log__.debug(*args, **kwargs)

    class LogShiftAdapter:
        def __init__(self, log, shift=8):
            self.log__ = log
            self.shift__ = shift
        def info(self, msg):
            self.log__.info(" " * self.shift__ + msg)

    log.info("Build info:")
    shlog = LogShiftAdapter(log)
    build_info.Log(shlog)

    log.info("\n")
    dt = datetime.utcnow()
    log.info("Current time (UTC ISO): " + dt.isoformat())
    log.info("\n")


#######################################################################################
#           Main
#######################################################################################
def main():
    binary_config = IgRepConConfig()
    log = CreateLogger()
    binary_config.CheckBinaries(log)
    parser, params = ParseCommandLineParams(log)
    if params.left_reads:
        CheckPairedReadsCorrectness(parser, params, log)
    else:
        CheckSingleReadsCorrectness(parser, params, log)
    CheckGeneralParamsCorrectness(parser, params, log)
    PrepareOutputDir(params)
    CreateFileLogger(params, log)
    PrintParams(params, log)
    PrintCommandLine(log)
    params.io = IgRepConIO(params.output, log)
    LogInfo(log)

    try:
        ig_phase_factory = PhaseFactory(params, log)
        ig_repertoire_constructor = PhaseManager(ig_phase_factory, params, log)
        ig_repertoire_constructor.Run()
        RemoveAuxFiles(params)
        PrintOutputFiles(params, log)
        log.info("\nThank you for using IgReC!")
    except (KeyboardInterrupt):
        log.info("\nIgReC was interrupted!")
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
