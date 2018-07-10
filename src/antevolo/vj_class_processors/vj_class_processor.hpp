#pragma once

#include "base_gene_class_processor.hpp"


namespace antevolo {
    class VJClassProcessor : public BaseGeneClassProcessor {
    public:
        VJClassProcessor(CloneSetWithFakesPtr clone_set_ptr,
                         const core::DecompositionClass& decomposition_class,
                         const AntEvoloConfig& config,
                         const AnnotatedCloneByReadConstructor& clone_by_read_constructor,
                         size_t current_fake_clone_index) :
                BaseGeneClassProcessor(clone_set_ptr,
                                       decomposition_class,
                                       config,
                                       clone_by_read_constructor,
                                       current_fake_clone_index) { }

        void CreateUniqueCDR3Map() override;
//        std::string WriteUniqueCDR3InFasta(core::DecompositionClass decomposition_class) override;
//        std::string GetGraphFname(core::DecompositionClass decomposition_class);

        EvolutionaryTree ProcessComponentWithKruskal(SparseGraphPtr hg_component, size_t component_id);
//        EvolutionaryTree ProcessComponentWithEdmonds(SparseGraphPtr hg_component, size_t component_id,
//                                                     const ShmModelEdgeWeightCalculator &edge_weight_calculator);
        std::vector<SparseGraphPtr> ComputeConnectedComponents() override;
    };
}