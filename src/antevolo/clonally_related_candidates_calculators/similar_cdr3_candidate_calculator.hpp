#pragma once

#include <decomposition.hpp>
#include <antevolo_config.hpp>
#include "candidate_calculator.hpp"
#include "../../graph_utils/sparse_graph.hpp"
#include <evolutionary_graph_utils/evolutionary_tree.hpp>
#include <evolutionary_graph_utils/evolutionary_edge_constructor.hpp>

namespace antevolo {
    class SimilarCDR3CandidateCalculator : public BaseCandidateCalculator {
        const AntEvoloConfig::AlgorithmParams &config_;
        const AntEvoloConfig::OutputParams &output_params_;
        size_t num_mismatches_;

        typedef std::map<std::string, std::vector<size_t>> UniqueCDR3IndexMap;

        UniqueCDR3IndexMap unique_cdr3s_map_;
        std::vector<std::string> unique_cdr3s_;

        SparseGraphPtr sparse_cdr_graph_;
        GraphComponentMap graph_component_;

        void Clear();

        std::string GetFastaFname(core::DecompositionClass decomposition_class);

    public:
        SimilarCDR3CandidateCalculator(const annotation_utils::CDRAnnotatedCloneSet &clone_set,
                                       const AntEvoloConfig::OutputParams &output_params,
                                       const AntEvoloConfig::AlgorithmParams &config) :
                BaseCandidateCalculator(clone_set),
                config_(config),
                output_params_(output_params),
                num_mismatches_(config.similar_cdr3s_params.num_mismatches) { }

        EvolutionaryEdgeConstructor* GetEdgeConstructor();
        void CreateUniqueCDR3Map(core::DecompositionClass decomposition_class);
        std::string WriteUniqueCDR3InFasta(core::DecompositionClass decomposition_class);
        std::string GetGraphFname(core::DecompositionClass decomposition_class);
        std::vector<SparseGraphPtr> ComputeCDR3HammingGraphs(std::string cdr_fasta, std::string cdr_graph);
        //ClonallyRelatedCandidates ComputeCandidatesForGraph(SparseGraphPtr hg_component, size_t component_id);
        EvolutionaryTree ComputeTreeForComponent(SparseGraphPtr hg_component, size_t component_id);
        void AddComponentToTheTree(SparseGraphPtr hg_component, size_t component_id, EvolutionaryTree& tree);
        //std::vector<ClonallyRelatedCandidates> ComputeCandidates(core::DecompositionClass decomposition_class);
        //std::vector<std::pair<size_t, size_t>> ComputeDirectedEdges(core::DecompositionClass decomposition_class);
    };
}