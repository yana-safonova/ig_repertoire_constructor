#include "antevolo_launch.hpp"

#include <read_archive.hpp>
#include <germline_utils/germline_db_generator.hpp>
#include <germline_db_labeler.hpp>
#include <vj_parallel_processor.hpp>
#include <read_labeler.hpp>
#include <cdr_output.hpp>
#include <evolutionary_graph_utils/evolutionary_tree_splitter.hpp>
#include <mutation_strategies/no_k_neighbours.hpp>
#include <evolutionary_graph_utils/one_child_fake_clones_filterer.hpp>
#include "antevolo_processor.hpp"

#include "antevolo_output_writer.hpp"
#include "annotated_clone_by_read_constructor.hpp"
#include "../fast_ig_tools/ig_trie_compressor.hpp"

#include "shm_model_utils/shm_model.hpp"

#include "shm_kmer_matrix_estimator/shm_kmer_matrix_estimator.hpp"
#include "posterior_distribution_calculator/posterior_distribution_calculator.hpp"
#include "kmer_matrix_exporter/kmer_matrix_exporter.hpp"
#include "stats_calculation/cluster_fillin_calculator.hpp"

#include "parallel_evolution/clonal_graph_writer.hpp"
//#include "parallel_evolution/mutation_map_writer.hpp"
//#include "parallel_evolution/trusted_shm_finder.hpp"
#include "parallel_evolution/clonal_graph_constructor.hpp"

#include "evolutionary_tree_annotation/annotated_tree_storage.hpp"

using namespace shm_kmer_matrix_estimator;

namespace antevolo {
    ShmModelEdgeWeightCalculator AntEvoloLaunch::ShmModelPosteriorCalculation(
            const annotation_utils::AnnotatedCloneSet<annotation_utils::AnnotatedClone>& clone_set)
    {
        ShmModel model(config_.input_params.shm_kmer_model_igh);

        ShmKmerMatrixEstimator statistics_estimator(config_.shm_config.mfp,
                                                    config_.shm_config.achp,
                                                    config_.shm_config.acrp);
        KmerMatrix fr_matrix, cdr_matrix;
        std::tie(fr_matrix, cdr_matrix) = statistics_estimator.calculate_mutation_statistics(clone_set);

        KmerMatrixExporter kmer_matrix_exporter;
        PosteriorDistributionCalculator posterior_distribution_calculator;
        ShmModel posterior_model(posterior_distribution_calculator.calculate(model, fr_matrix, cdr_matrix));

        AbstractMutationStrategyPtr
                mut_strategy(new NoKNeighboursMutationStrategy(config_.shm_config.mfp));
        ShmModelEdgeWeightCalculator edge_weight_calculator(posterior_model, std::move(mut_strategy));
        return edge_weight_calculator;
    }

    void AntEvoloLaunch::Launch() {
        INFO("AntEvolo starts");
//        INFO((config_.algorithm_params.model ? "edmonds" : "kruskal"));
//        INFO((config_.algorithm_params.compare ?
//              "comparing with respect to " +  config_.input_params.decomposition_rcm: "no comparing"));

        core::ReadArchive read_archive(config_.input_params.input_reads);
        if(config_.cdr_labeler_config.vj_finder_config.io_params.output_params.output_details.fix_spaces)
            read_archive.FixSpacesInHeaders();
        germline_utils::GermlineDbGenerator db_generator(config_.cdr_labeler_config.vj_finder_config.io_params.input_params.germline_input,
                                                         config_.cdr_labeler_config.vj_finder_config.algorithm_params.germline_params);
        INFO("Generation of DB for variable segments...");
        germline_utils::CustomGeneDatabase v_db = db_generator.GenerateVariableDb();
        INFO("Generation of DB for join segments...");
        germline_utils::CustomGeneDatabase j_db = db_generator.GenerateJoinDb();
        // todo: refactor code duplication
        INFO("CDR labeling for V gene segments");
        auto v_labeling = cdr_labeler::GermlineDbLabeler(v_db,
                                                         config_.cdr_labeler_config.cdrs_params).ComputeLabeling();
        INFO("CDR labeling for J gene segments");
        auto j_labeling = cdr_labeler::GermlineDbLabeler(j_db,
                                                         config_.cdr_labeler_config.cdrs_params).ComputeLabeling();
        INFO("Creation of labeled V and J databases");
        auto labeled_v_db = v_labeling.CreateFilteredDb();
        INFO("Labeled DB of V segments consists of " << labeled_v_db.size() << " records");
        auto labeled_j_db = j_labeling.CreateFilteredDb();
        INFO("Labeled DB of J segments consists of " << labeled_j_db.size() << " records");
        INFO("Alignment against VJ germline segments");
        vj_finder::VJParallelProcessor processor(read_archive,
                                                 config_.cdr_labeler_config.vj_finder_config.algorithm_params,
                                                 labeled_v_db, labeled_j_db,
                                                 config_.cdr_labeler_config.run_params.num_threads);
        vj_finder::VJAlignmentInfo alignment_info = processor.Process();
        INFO(alignment_info.NumVJHits() << " reads were aligned; " << alignment_info.NumFilteredReads() <<
             " reads were filtered out");
        if (config_.algorithm_params.compare && alignment_info.NumFilteredReads() > 0) {
            WARN("WARNING: Some reads were filtered out. EvoQuast mode assumes that all the reads have been cleaned before");
        }

        cdr_labeler::ReadCDRLabeler read_labeler(config_.cdr_labeler_config.shm_params, v_labeling, j_labeling);
        auto uncompressed_annotated_clone_set = read_labeler.CreateAnnotatedCloneSet(alignment_info);
        cdr_labeler::CDRLabelingWriter writer(config_.cdr_labeler_config.output_params,
                                              uncompressed_annotated_clone_set);

        //trie_compressor

        std::vector<seqan::Dna5String> clone_seqs;
        for (auto it = uncompressed_annotated_clone_set.cbegin(); it != uncompressed_annotated_clone_set.cend(); it++) {
            seqan::Dna5String clone_seq = it->Read().seq;
            clone_seqs.push_back(clone_seq);
        }
        INFO("Trie_compressor starts, " << uncompressed_annotated_clone_set.size() << " annotated sequences were created");
        auto indices = fast_ig_tools::Compressor::compressed_reads_indices(clone_seqs,
        fast_ig_tools::Compressor::Type::TrieCompressor);
        std::vector<size_t> abundances(uncompressed_annotated_clone_set.size());
        for (auto i : indices) {
            ++abundances[i];
        }
        annotation_utils::CDRAnnotatedCloneSet annotated_clone_set;
        size_t compressed_clone_index = 0;
        for (size_t i = 0; i < uncompressed_annotated_clone_set.size(); ++i) {
            if (abundances[i] != 0) {
                annotated_clone_set.AddClone(uncompressed_annotated_clone_set[i]);
                annotated_clone_set[compressed_clone_index].SetSize(abundances[i]);
                ++compressed_clone_index;
            }
        }
        INFO(annotated_clone_set.size() << " unique prefixes were created");

        //end trie_compressor

        writer.OutputCDRDetails();
        writer.OutputSHMs();

        AnnotatedCloneByReadConstructor clone_by_read_constructor(
                labeled_v_db,
                labeled_j_db,
                v_labeling,
                j_labeling,
                config_.cdr_labeler_config.vj_finder_config.algorithm_params,
                config_.cdr_labeler_config.shm_params);

        if (config_.algorithm_params.compare) {
            LaunchEvoQuast(annotated_clone_set);
        } else {
            LaunchDefault(clone_by_read_constructor,
                          annotated_clone_set,
                          read_archive.size());
        }

        //output_writer.OutputSHMForTrees();
        INFO("AntEvolo ends");
    }

    void AntEvoloLaunch::LaunchDefault(AnnotatedCloneByReadConstructor& clone_by_read_constructor,
                                       annotation_utils::CDRAnnotatedCloneSet& annotated_clone_set,
                                       size_t total_number_of_reads) {
        INFO("Tree construction starts");
        auto edge_weight_calculator = ShmModelPosteriorCalculation(annotated_clone_set);
        AntEvoloProcessor antevolo_processor = AntEvoloProcessor(config_,
                                                                 annotated_clone_set,
                                                                 clone_by_read_constructor,
                                                                 total_number_of_reads,
                                                                 edge_weight_calculator);

        auto tree_storage = antevolo_processor.ConstructClonalTreesHG();
        //auto tree_storage = antevolo_processor.GetCDR3Stats();

        auto final_clone_set = antevolo_processor.GetCloneSetWithFakes();
        INFO("Evolutionary directions for " << tree_storage.size() << " clonal lineages were created");
        INFO("Computation of evolutionary statistics");
        // todo: add refactoring!!!
        EvolutionaryTreeStorage connected_tree_storage;
        OneChildFakeClonesFilterer fakes_filterer(config_.algorithm_params.edge_construction_params);
        for(auto it = tree_storage.cbegin(); it != tree_storage.cend(); it++) {
            ConnectedTreeSplitter tree_splitter;
            auto connected_trees = tree_splitter.Split(*it);
            for(auto it2 = connected_trees.begin(); it2!= connected_trees.end(); it2++) {
                auto filtered_tree = fakes_filterer.FilterOneChildFakes(*it2);
                connected_tree_storage.Add(filtered_tree);
            }
        }
        INFO(tree_storage.size() << " clonal lineages were splitted into " << connected_tree_storage.size() <<
                                 " connected trees");

        AnalyzeParallelEvolution(connected_tree_storage);

        AnnotatedTreeStorage annotated_storage;
        for(auto it = connected_tree_storage.cbegin(); it != connected_tree_storage.cend(); it++) {
            annotated_storage.AddAnnotatedTree(*it);
        }
        INFO("Annotation for " << annotated_storage.size() << " clonal trees was computed");

        AntEvoloOutputWriter output_writer(config_.output_params, annotated_storage);
        output_writer.OutputTreeStats();
        output_writer.OutputSHMForTrees();

        output_writer.OutputCleanedSequences(final_clone_set);
        INFO("Cleaned sequences were written to " << config_.output_params.output_dir << "/cleaned_sequences.fa");

        //for (auto it = connected_tree_storage.cbegin(); it != connected_tree_storage.cend(); it++) {
        for (auto it = tree_storage.cbegin(); it != tree_storage.cend(); it++) {
            output_writer.WriteTreeInFile(config_.output_params.tree_dir, *it);
            output_writer.WriteTreeVerticesInFile(config_.output_params.vertex_dir, *it);
            //TRACE(i + 1 << "-th clonal tree was written to " << tree.Get);
        }
        output_writer.WriteRcmFromStorageInFile(config_.output_params.output_dir, connected_tree_storage);

        INFO("Clonal trees were written to " << config_.output_params.tree_dir);
    };

    void AntEvoloLaunch::AnalyzeParallelEvolution(const EvolutionaryTreeStorage& trees) {
        if(!config_.algorithm_params.parallel_evolution_params.enable_parallel_shms_finder)
            return;
        INFO("Parallel evolution finder starts");
        for(auto it = trees.cbegin(); it != trees.cend(); it++) {
            ClonalGraph cgraph = ClonalGraphConstructor(*it).ComputeClonalGraph();
            UniqueSHMCalculator shm_calc(it->GetCloneSet(), *it);
            for(auto e = cgraph.cbegin(); e != cgraph.cend(); e++) {
                auto src_nodes = e->second;
                for(auto src = src_nodes.begin(); src != src_nodes.end(); src++)
                    shm_calc.AddSHMsFromEdge(e->first, *src);
            }
            TreeSHMMap shm_map = shm_calc.GetSHMMap();
            ClonalGraphWriter graph_writer(cgraph);
            graph_writer(path::append_path(config_.output_params.parallel_shm_output.parallel_bulges_dir,
                                           it->GetTreeOutputFname("") + ".dot"));
        }
    }

    void AntEvoloLaunch::LaunchEvoQuast(const annotation_utils::CDRAnnotatedCloneSet& clone_set) {
        INFO("EvoQuast starts");
        boost::unordered_map<std::string, size_t> read_name_to_index;
        for (size_t i = 0; i < clone_set.size(); ++i) {
            read_name_to_index[clone_set[i].Read().name] = i;
        }

        std::vector<boost::unordered_set<size_t>> clusters = ReadClusters(read_name_to_index);

        auto fillin_calculator = ClusterFillinCalculator(clone_set, config_);
        std::vector<size_t> clusters_edge_numbers(clusters.size());

        size_t fracturing_stat = clusters.size();
        std::vector<size_t> clusters_Es(clusters.size());
#pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < clusters.size(); ++i) {
            clusters_Es[i] = fillin_calculator.CalculateE(clusters[i]);
        }

        for (size_t i = 0; i < clusters.size(); ++i) {
            fracturing_stat += 2 * clusters_Es[i];
        }


        INFO("Fracturing statistic: " << fracturing_stat);
        INFO("EvoQuast ends");
    };

    std::vector<boost::unordered_set<size_t>> AntEvoloLaunch::ReadClusters(
            const boost::unordered_map<std::string, size_t>& read_name_to_index) {
        std::vector<boost::unordered_set<size_t>> clusters;
        std::cout << config_.input_params.decomposition_rcm << std::endl;
        std::fstream in(config_.input_params.decomposition_rcm);
        VERIFY_MSG(in.is_open(), std::string("RCM file is not opened!"));
        std::string line;
        while (getline(in, line) && !line.empty()) {
            size_t tab = line.find('\t');
            std::string read_name = line.substr(0, tab);
            std::stringstream ss;
            ss << line.substr(tab+1, line.length());
            size_t cluster_num;
            ss >> cluster_num;
            while (clusters.size() < cluster_num + 1) {
                clusters.push_back(boost::unordered_set<size_t>());
            }
            clusters[cluster_num].insert(read_name_to_index.find(read_name)->second);
        }
        in.close();
        return clusters;
    }
}
