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
        ShmModel& model_;
        const AnnotatedCloneByReadConstructor& clone_by_read_constructor_;

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
        VJClassProcessor(const annotation_utils::CDRAnnotatedCloneSet& clone_set,
                         const AntEvoloConfig::OutputParams &output_params,
                         const AntEvoloConfig::AlgorithmParams &config,
                         ShmModel& model,
                         const AnnotatedCloneByReadConstructor& clone_by_read_constructor) :
                BaseCandidateCalculator(clone_set),
                config_(config),
                output_params_(output_params),
                num_mismatches_(config.similar_cdr3s_params.num_mismatches),
                model_(model),
                clone_by_read_constructor_(clone_by_read_constructor) { }

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