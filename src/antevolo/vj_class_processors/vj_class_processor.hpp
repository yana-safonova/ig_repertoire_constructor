#pragma once

#include <decomposition.hpp>
#include <antevolo_config.hpp>
#include "base_vj_class_processor.hpp"
#include "cdr3_hamming_graph_connected_components_processors/kruskal_cdr3_hg_cc_processor.hpp"

namespace antevolo {
    class VJClassProcessor : public BaseCandidateCalculator {

        const AntEvoloConfig::AlgorithmParams &config_;
        const AntEvoloConfig::OutputParams &output_params_;
        size_t num_mismatches_;
        const AnnotatedCloneByReadConstructor& clone_by_read_constructor_;
        size_t& current_fake_clone_index_;
        size_t& reconstructed_;
        size_t& rejected_;

        typedef std::map<std::string, std::vector<size_t>> UniqueCDR3IndexMap;
        typedef std::map<std::string, size_t> CDR3ToIndexMap;
        typedef boost::associative_property_map<std::map<size_t, size_t>> AP_map;

        UniqueCDR3IndexMap unique_cdr3s_map_;
        CDR3ToIndexMap cdr3_to_old_index_map_;
        std::vector<std::string> unique_cdr3s_;

        SparseGraphPtr sparse_cdr_graph_;
        GraphComponentMap graph_component_map_;



        void Clear();

        std::string GetFastaFname(core::DecompositionClass decomposition_class);

    public:
        VJClassProcessor(CloneSetWithFakesPtr clone_set,
                         const AntEvoloConfig::OutputParams &output_params,
                         const AntEvoloConfig::AlgorithmParams &config,
                         const AnnotatedCloneByReadConstructor& clone_by_read_constructor,
                         size_t& current_fake_clone_index,
                         size_t& reconstructed,
                         size_t& rejected) :
                BaseCandidateCalculator(clone_set),
                config_(config),
                output_params_(output_params),
                num_mismatches_(config.similar_cdr3s_params.num_mismatches),
                clone_by_read_constructor_(clone_by_read_constructor),
                current_fake_clone_index_(current_fake_clone_index),
                reconstructed_(reconstructed),
                rejected_(rejected) { }

        //EvolutionaryEdgeConstructor* GetEdgeConstructor();
        void CreateUniqueCDR3Map(core::DecompositionClass decomposition_class);
        std::string WriteUniqueCDR3InFasta(core::DecompositionClass decomposition_class);
        std::string GetGraphFname(core::DecompositionClass decomposition_class);
        std::vector<SparseGraphPtr> ComputeCDR3HammingGraphs(std::string cdr_fasta, std::string cdr_graph);
        //void AddUndirectedForestToTheTree(SparseGraphPtr hg_component, size_t component_id, EvolutionaryTree& tree,
        //                                  boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges);
        //void AddComponentToTheTree(SparseGraphPtr hg_component, size_t component_id, EvolutionaryTree& tree);
        //void SetUndirectedComponentsParentEdges(SparseGraphPtr hg_component,
        //                                        size_t component_id,
        //                                        EvolutionaryTree& tree,
        //                                        boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges);
        //void SetDirections(boost::unordered_set<size_t> vertices_nums,
        //                                        EvolutionaryTree& tree,
        //                                        boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges);
        EvolutionaryTree AddComponent(SparseGraphPtr hg_component, size_t component_id);
    };
}