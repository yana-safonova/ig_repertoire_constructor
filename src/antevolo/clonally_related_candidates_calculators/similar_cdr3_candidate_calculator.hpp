#pragma once

#include <decomposition.hpp>
#include <antevolo_config.hpp>
#include "candidate_calculator.hpp"
#include "../../graph_utils/sparse_graph.hpp"

namespace antevolo {
    class SimilarCDR3CandidateCalculator : public BaseCandidateCalculator {
        const AntEvoloConfig::OutputParams &output_params_;
        size_t num_mismatches_;

        typedef std::map<std::string, std::vector<size_t>> UniqueCDR3IndexMap;

        UniqueCDR3IndexMap unique_cdr3s_map_;
        std::vector<std::string> unique_cdr3s_;

        SparseGraphPtr sparse_cdr_graph_;
        GraphComponentMap graph_component_;

        void Clear();

        void CreateUniqueCDR3Map(core::DecompositionClass decomposition_class);

        std::string GetFastaFname(core::DecompositionClass decomposition_class);

        std::string WriteUniqueCDR3InFasta(core::DecompositionClass decomposition_class);

        std::string GetGraphFname(core::DecompositionClass decomposition_class);

        std::vector<SparseGraphPtr> ComputeCDR3HammingGraphs(std::string cdr_fasta, std::string cdr_graph);

        ClonallyRelatedCandidates ComputeCandidatesForGraph(SparseGraphPtr hg_component, size_t component_id);

    public:
        SimilarCDR3CandidateCalculator(const annotation_utils::CDRAnnotatedCloneSet &clone_set,
                                       const AntEvoloConfig::OutputParams &output_params,
                                       size_t num_mismatches) :
                BaseCandidateCalculator(clone_set),
                output_params_(output_params),
                num_mismatches_(num_mismatches) { }

        std::vector<ClonallyRelatedCandidates> ComputeCandidates(core::DecompositionClass decomposition_class);
    };
}