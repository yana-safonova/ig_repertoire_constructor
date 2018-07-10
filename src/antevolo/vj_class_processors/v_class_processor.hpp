#pragma once

#include "base_gene_class_processor.hpp"
#include <cdr3_hamming_graph_component_info.hpp>

namespace antevolo {
    class VClassProcessor : public BaseGeneClassProcessor {
        std::vector<size_t> jdifference_positions_;

        EvolutionaryTree ProcessComponentWithEdmonds(SparseGraphPtr hg_component, size_t component_id,
                                                     const ShmModelEdgeWeightCalculator &edge_weight_calculator);
    public:

        VClassProcessor(CloneSetWithFakesPtr clone_set_ptr,
                        const core::DecompositionClass& decomposition_class,
                        const AntEvoloConfig& config,
                        const AnnotatedCloneByReadConstructor& clone_by_read_constructor,
                        size_t current_fake_clone_index) :
                BaseGeneClassProcessor(clone_set_ptr,
                                       decomposition_class,
                                       config,
                                       clone_by_read_constructor,
                                       current_fake_clone_index) {
            auto chain = clone_set_ptr->operator[](*decomposition_class_.cbegin()).ChainType().Chain();
            if (chain == germline_utils::ImmuneChainType::HeavyIgChain) {
                jdifference_positions_ = {17, 18, 19, 22, 25, 26, 27}; // FIXME: move to config!!
            } else if (chain == germline_utils::ImmuneChainType::KappaIgChain) {
                jdifference_positions_ = {12, 13, 14, 15, 16, 22, 23, 24};
            } else if (chain == germline_utils::ImmuneChainType::LambdaIgChain) {
                jdifference_positions_ = {12, 15, 19, 22, 23, 24, 25};
            } else {
                jdifference_positions_ = {};
            }
            VERIFY(jdifference_positions_.size() != 0);
        }

        void CreateUniqueCDR3Map() override;
        std::vector<SparseGraphPtr> ComputeConnectedComponents() override;

        void ChangeJgene(const germline_utils::ImmuneGene &v_gene,
                         const germline_utils::ImmuneGene &j_gene);
        void ChangeJgeneToMax(CDR3HammingGraphComponentInfo hamming_graph_info);
    };
}