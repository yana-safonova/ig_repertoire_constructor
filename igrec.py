#!/usr/bin/env python2

import os
import shutil
import sys
import logging

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
class PhaseNames:
    def __init__(self):
        self.__pair_reads_merger = 'pair_reads_merger'
        self.__vj_alignment = 'vj_alignment'
        self.__trie_compressor = 'trie_compressor'
        self.__graph_construction = 'graph_constructor'
        self.__dsf = 'dsf'
        self.__consensus_constructor = 'consensus_constructor'
        self.__compress_equal_clusters = 'compress_equal_clusters'
        self.__remove_low_abundance_reads = 'remove_low_abundance_reads'
        self.__phase_order = [self.__pair_reads_merger,
                              self.__vj_alignment,
                              self.__trie_compressor,
                              self.__graph_construction,
                              self.__dsf,
                              self.__consensus_constructor,
                              self.__compress_equal_clusters,
                              self.__remove_low_abundance_reads]
        self.__long_names = {'pair_reads_merger': 'Pair reads merging',
                             'vj_alignment' : 'VJ Alignment',
                             'trie_compressor' : 'Trie Compressor',
                             'graph_constructor' : 'Graph Constructor',
                             'dsf' : 'Dense Subgraph Finder',
                             'consensus_constructor' : 'Consensus Constructor',
                             'compress_equal_clusters' : 'Compress Equal Final Clusters',
                             'remove_low_abundance_reads' : 'Low Abundant Clusters Remover'}

    def __iter__(self):
        for sname in self.__phase_order:
            yield sname

    def __len__(self):
        return len(self.__phase_order)

    def GetPhaseNameBy(self, index):
        return self.__phase_order[index]

    def GetPhaseIndex(self, phase_name):
        for i in range(len(self)):
            if self.GetPhaseNameBy(i) == phase_name:
                return i
        return -1

    def PhaseIsPairReadsMerger(self, phase_name):
        return phase_name == self.__pair_reads_merger

    def GetPairReadMergerLongName(self):
        return self.__long_names[self.__pair_reads_merger]

    def PhaseIsVJAlignment(self, phase_name):
        return phase_name == self.__vj_alignment

    def GetVJAlignmentLongName(self):
        return self.__long_names[self.__vj_alignment]

    def PhaseIsTrieCompressor(self, phase_name):
        return phase_name == self.__trie_compressor

    def GetTrieCompressorLongName(self):
        return self.__long_names[self.__trie_compressor]

    def PhaseIsGraphConstructor(self, phase_name):
        return phase_name == self.__graph_construction

    def GetGraphConstructionLongName(self):
        return self.__long_names[self.__graph_construction]

    def PhaseIsDSF(self, phase_name):
        return phase_name == self.__dsf

    def GetDSFLongName(self):
        return self.__long_names[self.__dsf]

    def PhaseIsConsensusConstructor(self, phase_name):
        return phase_name == self.__consensus_constructor

    def GetConsensusConstructorLongName(self):
        return self.__long_names[self.__consensus_constructor]

    def PhaseIsCompressEqualClusters(self, phase_name):
        return phase_name == self.__compress_equal_clusters

    def GetCompressEqualClustersName(self):
        return self.__long_names[self.__compress_equal_clusters]

    def PhaseIsRemoveLowAbundanceReads(self, phase_name):
        return phase_name == self.__remove_low_abundance_reads

    def GetRemoveLowAbundanceReadsName(self):
        return self.__long_names[self.__remove_low_abundance_reads]

###########
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
        self.path_to_dsf = os.path.join(home_directory, 'build/release/bin/dense_sgraph_finder')
        self.path_to_germline = os.path.join(home_directory, "data/germline")

    def __init__(self):
        self.__initBinaryPaths()

    def CheckBinaries(self, log):
        phase_names = PhaseNames()
        if not os.path.exists(self.path_to_pair_reads_merger):
            log.info("ERROR: Binary file of " + phase_names.GetPairReadMergerLongName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_vj_aligner):
            log.info("ERROR: Binary file of " + phase_names.GetVJAlignmentLongName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_trie_compressor):
            log.info("ERROR: Binary file of " + phase_names.GetTrieCompressorLongName() + " (" + self.path_to_trie_compressor +") was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.run_report_supernodes):
            log.info("ERROR: Binary file of " + phase_names.GetTrieCompressorLongName() +  " (" + self.run_report_supernodes + ") was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.run_triecmp_to_repertoire):
            log.info("ERROR: Binary file of " + phase_names.GetTrieCompressorLongName() + " (" + self.run_triecmp_to_repertoire + ") was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_graph_constructor):
            log.info("ERROR: Binary file of " + phase_names.GetGraphConstructionLongName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_dsf):
            log.info("ERROR: Binary file of " + phase_names.GetDSFLongName() + " was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.path_to_consensus_constructor):
            log.info("ERROR: Binary file of " + phase_names.GetConsensusConstructorLongName() + " (" + self.path_to_consensus_constructor + ") was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.run_rcm_recoverer):
            log.info("ERROR: Binary file of " + phase_names.GetConsensusConstructorLongName() + " (" + self.run_rcm_recoverer + ") was not found\n")
            ErrorMessagePrepareCfg(log)
            sys.exit(1)
        if not os.path.exists(self.run_compress_equal_clusters):
            log.info("ERROR: Binary file of " + phase_names.GetCompressEqualClustersName() + " was not found\n")
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

    def __init__(self, output_dir, log):
        self.__log = log
        self.__initVJFinderOutput(output_dir)
        self.__initCompressorOutput(output_dir)
        self.sw_graph = os.path.join(output_dir, "sw.graph")
        self.__initDSFOutput(output_dir)
        self.__initFinalOutput(output_dir)
        self.final_stripped_clusters_fa = os.path.join(output_dir, 'final_repertoire_large.fa')
        self.__initCompressEqualClusters(output_dir)

    def CheckCroppedReadsExistance(self):
        if not os.path.exists(self.cropped_reads):
            self.__log.info("ERROR: File containing cleaned Ig-Seq reads was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckBadReadsExistance(self):
        if not os.path.exists(self.bad_reads):
            self.__log.info("ERROR: File containing contaminated reads (not Ig-Seq) was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckVJAlignmentInfoExistance(self):
        if not os.path.exists(self.vj_alignment_info):
            self.__log.info("ERROR: File containing VJ alignment info was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckCompressedReadsExistance(self):
        if not os.path.exists(self.compressed_reads):
            self.__log.info("ERROR: File containing compressed reads was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckCroppedCompressedMapExistance(self):
        if not os.path.exists(self.map_file):
            self.__log.info("ERROR: File containing map from cleaned reads to compressed reads was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckSupernodesExistance(self):
        if not os.path.exists(self.supernodes_file):
            self.__log.info("ERROR: File containing super-reads was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckSWGraphExistance(self):
        if not os.path.exists(self.sw_graph):
            self.__log.info("ERROR: File containing Smith-Waterman graph was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckDenseSubgraphDecompositionExistance(self):
        if not os.path.exists(self.dense_sgraph_decomposition):
            self.__log("ERROR: File containing dense subgraph decomposition was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckUncompressedFinalClustersExistance(self):
        if not os.path.exists(self.uncompressed_final_clusters_fa):
            self.__log("ERROR: File containing uncompressed antibody clusters of final repertoire was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckCompressedFinalClustersExistance(self):
        if not os.path.exists(self.compressed_final_clusters_fa):
            self.__log("ERROR: File containing compressed antibody clusters of final repertoire was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckUncompressedFinalRCMExistance(self):
        if not os.path.exists(self.uncompressed_final_rcm):
            self.__log("ERROR: File containing RCM of uncompressed final repertoire was not found")
            SupportInfo(self.__log)
            sys.exit(1)
    def CheckCompressedFinalRCMExistance(self):
        if not os.path.exists(self.compressed_final_rcm):
            self.__log("ERROR: File containing RCM of compressed final repertoire was not found")
            SupportInfo(self.__log)
            sys.exit(1)

    def CheckFinalStrippedClustersExistance(self):
        if not os.path.exists(self.final_stripped_clusters_fa):
            self.__log("ERROR: File containing large antibody clusters of final repertoire was not found")
            SupportInfo(self.__log)
            sys.exit(1)

#######################################################################################
#           Phases
#######################################################################################
class Phase:
    def __init__(self, long_name, log):
        self._long_name = long_name
        self._log = log

    def PrintStartMessage(self):
        self._log.info("==== " + self._long_name + " starts\n")

    def Run(self):
        print "This method should be overloaded"

    def PrintFinishMessage(self):
        self._log.info("\n==== " + self._long_name + " finished")

###########

class PairReadMerger(Phase):
    def __init__(self, params, log):
        Phase.__init__(self, PhaseNames().GetPairReadMergerLongName(), log)
        self.__params = params

    def __CheckInputExistance(self):
        if not os.path.exists(self.__params.left_reads):
            self._log.info("ERROR: Input left reads " + self.__params.left_reads + " were not found")
            SupportInfo(self._log)
            sys.exit(1)
        if not os.path.exists(self.__params.right_reads):
            self._log.info("ERROR: Input right reads " + self.__params.right_reads + " were not found")
            SupportInfo(self._log)
            sys.exit(1)

    def __CheckOutputExistance(self):
        if not os.path.exists(self.__params.single_reads):
            self._log.info("ERROR: Input reads " + self.__params.single_reads + " were not found")
            SupportInfo(self._log)
            sys.exit(1)

    def Run(self):
        self.__CheckInputExistance()
        command_line = "%s %s %s %s" % (IgRepConConfig().run_pair_reads_merger,
                                        self.__params.left_reads,
                                        self.__params.right_reads,
                                        self.__params.single_reads)
        support.sys_call(command_line, self._log)

    def PrintOutputFiles(self):
        self.__CheckOutputExistance()
        self._log.info("\nOutput files: ")
        self._log.info("  * Merged reads were written to " + self.__params.single_reads)

###########
class VJAlignmentPhase(Phase):
    def __init__(self, params, log):
        Phase.__init__(self, PhaseNames().GetVJAlignmentLongName(), log)
        self.__params = params

    def __CheckInputExistance(self):
        if not os.path.exists(self.__params.single_reads):
            self._log.info("ERROR: Input reads " + self.__params.single_reads + " were not found")
            SupportInfo(self._log)
            sys.exit(1)

    def __CheckOutputExistance(self):
        self.__params.io.CheckCroppedReadsExistance()
        if not self.__params.no_alignment:
            self.__params.io.CheckBadReadsExistance()
            self.__params.io.CheckVJAlignmentInfoExistance()

    def Run(self):
        self.__CheckInputExistance()
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
            cwd = os.getcwd()
            os.chdir(home_directory)
            support.sys_call(command_line, self._log)
            os.chdir(cwd)
        else:
            self._log.info("VJ Finder stage skipped")
            self.__params.io.cropped_reads = self.__params.single_reads

    def PrintOutputFiles(self):
        self.__CheckOutputExistance()
        if not self.__params.no_alignment:
            self._log.info("\nOutput files: ")
            self._log.info("  * Cleaned Ig-Seq reads were written to " + self.__params.io.cropped_reads)
            self._log.info("  * Contaminated (not Ig-Seq) reads were written to " + self.__params.io.bad_reads)
            self._log.info("  * VJ alignment output was written to " + self.__params.io.vj_alignment_info)

###########
class TrieCompressionPhase(Phase):
    def __init__(self, params, log):
        Phase.__init__(self, PhaseNames().GetTrieCompressorLongName(), log)
        self.__params = params

    def __CheckInputExistance(self):
        self.__params.io.CheckCroppedReadsExistance()

    def __CheckOutputExistance(self):
        self.__params.io.CheckCompressedReadsExistance()
        self.__params.io.CheckCroppedCompressedMapExistance()
        self.__params.io.CheckSupernodesExistance()

    def Run(self):
        self.__CheckInputExistance()
        command_line = IgRepConConfig().run_trie_compressor + " -i " + self.__params.io.cropped_reads + \
                    " -o " + self.__params.io.compressed_reads + " -m " + self.__params.io.map_file + " -Toff"
        support.sys_call(command_line, self._log)

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


    def PrintOutputFiles(self):
        self.__CheckOutputExistance()
        self._log.info("\nOutput files:")
        self._log.info("  * Compressed reads were written to " + self.__params.io.compressed_reads)
        self._log.info("  * Super reads were written to " + self.__params.io.supernodes_file)

###########
class GraphConstructionPhase(Phase):
    def __init__(self, params, log):
        Phase.__init__(self, PhaseNames().GetGraphConstructionLongName(), log)
        self.__params = params

    def __CheckInputExistance(self):
        self.__params.io.CheckCompressedReadsExistance()

    def __CheckOutputExistance(self):
        self.__params.io.CheckSWGraphExistance()

    def Run(self):
        self.__CheckInputExistance()
        command_line = IgRepConConfig().run_graph_constructor + " -i " + self.__params.io.compressed_reads + \
                       " -o " + self.__params.io.sw_graph + " -t " + str(self.__params.num_threads) + \
                       " --tau=" + str(self.__params.max_mismatches) + " -A" + " -Toff"
        support.sys_call(command_line, self._log)

    def PrintOutputFiles(self):
        self.__CheckOutputExistance()
        self._log.info("\nOutput files:")
        self._log.info("  * Smith-Waterman graph was written to " + self.__params.io.sw_graph)

###########
class DSFPhase(Phase):
    def __init__(self, params, log):
        Phase.__init__(self, PhaseNames().GetDSFLongName(), log)
        self.__params = params

    def __CheckInputExistance(self):
        self.__params.io.CheckSWGraphExistance()

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

    def __CheckOutputExistance(self):
        self.__params.io.CheckDenseSubgraphDecompositionExistance()

    def Run(self):
        self.__CheckInputExistance()
        dense_subgraph_finder.main(self.__GetDSFParams(), self.__params.log_filename)

    def PrintOutputFiles(self):
        self.__CheckOutputExistance()
        self._log.info("\nOutput files:")
        self._log.info("  * Dense subgraph decomposition was written to " +
                       self.__params.io.dense_sgraph_decomposition)

###########
class ConsensusConstructionPhase(Phase):
    def __init__(self, params, log):
        Phase.__init__(self, PhaseNames().GetConsensusConstructorLongName(), log)
        self.__params = params

    def __CheckInputExistance(self):
        self.__params.io.CheckCompressedReadsExistance()
        self.__params.io.CheckDenseSubgraphDecompositionExistance()
        self.__params.io.CheckCroppedReadsExistance()
        self.__params.io.CheckCroppedCompressedMapExistance()

    def __CheckOutputExistance(self):
        self.__params.io.CheckUncompressedFinalClustersExistance()
        self.__params.io.CheckUncompressedFinalRCMExistance()

    def Run(self):
        self.__CheckInputExistance()
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
        support.sys_call(command_line, self._log)


    def PrintOutputFiles(self):
        self.__CheckOutputExistance()
        self._log.info("\nOutput files:")
        self._log.info("  * Antibody clusters of uncompressed final repertoire were written to " +
                       self.__params.io.uncompressed_final_clusters_fa)
        self._log.info("  * Read-cluster map of uncompressed final repertoire was written to " +
                       self.__params.io.uncompressed_final_rcm)

class CompressEqualClusters(Phase):
    def __init__(self, params, log):
        Phase.__init__(self, PhaseNames().GetCompressEqualClustersName(), log)
        self.__params = params

    def __CheckInputExistance(self):
        self.__params.io.CheckUncompressedFinalClustersExistance()
        self.__params.io.CheckUncompressedFinalRCMExistance()

    def __CheckOutputExistance(self):
        self.__params.io.CheckCompressedFinalClustersExistance()
        self.__params.io.CheckCompressedFinalRCMExistance()

    def Run(self):
        self.__CheckInputExistance()
        command_line = "%s %s %s -T %s -m %s -r %s -R %s" % (IgRepConConfig().run_compress_equal_clusters,
                                                             self.__params.io.uncompressed_final_clusters_fa,
                                                             self.__params.io.compressed_final_clusters_fa,
                                                             self.__params.io.tmp_compressed_clusters_fa,
                                                             self.__params.io.tmp_compressed_clusters_map,
                                                             self.__params.io.uncompressed_final_rcm,
                                                             self.__params.io.compressed_final_rcm)
        support.sys_call(command_line, self._log)


    def PrintOutputFiles(self):
        self.__CheckOutputExistance()
        self._log.info("\nOutput files:")
        self._log.info("  * Equal output clusters joined " +
                       self.__params.io.compressed_final_clusters_fa)

class RemoveLowAbundanceReadsPhase(Phase):
    def __init__(self, params, log):
        Phase.__init__(self, PhaseNames().GetRemoveLowAbundanceReadsName(), log)
        self.__params = params

    def __CheckInputExistance(self):
        self.__params.io.CheckCompressedFinalClustersExistance()

    def __CheckOutputExistance(self):
        self.__params.io.CheckFinalStrippedClustersExistance()

    def Run(self):
        self.__CheckInputExistance()
        command_line = "%s %s %s --limit=%d" % (IgRepConConfig().run_report_supernodes,
                                                self.__params.io.compressed_final_clusters_fa,
                                                self.__params.io.final_stripped_clusters_fa,
                                                self.__params.min_cluster_size)
        support.sys_call(command_line, self._log)


    def PrintOutputFiles(self):
        self.__CheckOutputExistance()
        self._log.info("\nOutput files:")
        self._log.info("  * Highly abundant antibody clusters of final repertoire were written to " +
                       self.__params.io.final_stripped_clusters_fa)

###########
class PhaseFactory:
    def __init__(self, phase_names, params, log):
        self.__phase_names = phase_names
        self.__entry_point = params.entry_point
        self.__params = params
        self.__log = log

    def __CreatePhaseByName(self, phase_name):
        if self.__phase_names.PhaseIsPairReadsMerger(phase_name):
            return PairReadMerger(self.__params, self.__log)
        elif self.__phase_names.PhaseIsVJAlignment(phase_name):
            return VJAlignmentPhase(self.__params, self.__log)
        elif self.__phase_names.PhaseIsTrieCompressor(phase_name):
            return TrieCompressionPhase(self.__params, self.__log)
        elif self.__phase_names.PhaseIsGraphConstructor(phase_name):
            return GraphConstructionPhase(self.__params, self.__log)
        elif self.__phase_names.PhaseIsDSF(phase_name):
            return DSFPhase(self.__params, self.__log)
        elif self.__phase_names.PhaseIsConsensusConstructor(phase_name):
            return ConsensusConstructionPhase(self.__params, self.__log)
        elif self.__phase_names.PhaseIsCompressEqualClusters(phase_name):
            return CompressEqualClusters(self.__params, self.__log)
        elif self.__phase_names.PhaseIsRemoveLowAbundanceReads(phase_name):
            return RemoveLowAbundanceReadsPhase(self.__params, self.__log)

    def CreatePhases(self):
        phase_list = list()
        first_phase_index = self.__phase_names.GetPhaseIndex(self.__entry_point)
        if first_phase_index == -1:
            self.__log.info("Incorrect name of entry-point")
            sys.exit(1)
        for i in range(first_phase_index, len(self.__phase_names )):
            phase_list.append(self.__CreatePhaseByName(self.__phase_names.GetPhaseNameBy(i)))
        return phase_list

############
class PhaseManager:
    def __init__(self, phase_factory, params, log):
        self.__params = params
        self.__log = log
        self.__phase_factory = phase_factory
        self.__phases = self.__phase_factory.CreatePhases()

    def __RunSinglePhase(self, phase_index):
            self.__phases[phase_index].PrintStartMessage()
            self.__phases[phase_index].Run()
            self.__phases[phase_index].PrintOutputFiles()
            self.__phases[phase_index].PrintFinishMessage()

    def __PrintPhaseDelimeter(self):
        self.__log.info("\n============================================\n")

    def Run(self, start_phase=0):
        self.__RunSinglePhase(start_phase)
        for i in range(start_phase + 1, len(self.__phases) - 1):
            self.__PrintPhaseDelimeter()
            self.__RunSinglePhase(i)
        if len(self.__phases) - start_phase != 1:
            self.__PrintPhaseDelimeter()
            self.__RunSinglePhase(len(self.__phases) - 1)

#######################################################################################
#           IO routines
#######################################################################################

def CreateLogger():
    log = logging.getLogger('igrec')
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
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
            setattr(namespace, "loci", "all")
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
                          default=PhaseNames().GetPhaseNameBy(0),
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
    log.info("  Entry point:\t\t\t" + params.entry_point)

def CreateFileLogger(params, log):
    params.log_filename = os.path.join(params.output, "igrec.log")
    if os.path.exists(params.log_filename):
        os.remove(params.log_filename)
    log_handler = logging.FileHandler(params.log_filename, mode='a')
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

    try:
        ig_phase_factory = PhaseFactory(PhaseNames(), params, log)
        ig_repertoire_constructor = PhaseManager(ig_phase_factory, params, log)
        if params.left_reads:
            ig_repertoire_constructor.Run(start_phase=0)
        else:
            ig_repertoire_constructor.Run(start_phase=1)
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
