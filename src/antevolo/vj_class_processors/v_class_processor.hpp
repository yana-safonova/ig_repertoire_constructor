#pragma once

#include "base_gene_class_processor.hpp"


namespace antevolo {
    class VClassProcessor : public BaseGeneClassProcessor {



    public:
        VClassProcessor(CloneSetWithFakesPtr clone_set_ptr,
                         const AntEvoloConfig& config,
                         const AnnotatedCloneByReadConstructor& clone_by_read_constructor,
                         size_t current_fake_clone_index) :
                BaseGeneClassProcessor(clone_set_ptr,
                                       config,
                                       clone_by_read_constructor,
                                       current_fake_clone_index) { }

//        void CreateUniqueCDR3Map(core::DecompositionClass decomposition_class);
        std::string WriteUniqueCDR3InFasta(core::DecompositionClass decomposition_class) override;
//        std::string GetGraphFname(core::DecompositionClass decomposition_class);

        EvolutionaryTree ProcessComponentWithKruskal(SparseGraphPtr hg_component, size_t component_id);
//        EvolutionaryTree ProcessComponentWithEdmonds(SparseGraphPtr hg_component, size_t component_id,
//                                                     const ShmModelEdgeWeightCalculator &edge_weight_calculator);
        vector<SparseGraphPtr> ComputeConnectedComponents(const core::DecompositionClass& decomposition_class) override;
    };
}