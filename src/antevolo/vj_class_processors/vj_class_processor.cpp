#include "vj_class_processor.hpp"
#include <cdr3_hamming_graph_connected_components_processors/edmonds_cdr3_hg_cc_processor.hpp>


namespace antevolo {

    EvolutionaryTree VJClassProcessor::ProcessComponentWithKruskal(SparseGraphPtr hg_component, size_t component_id) {

        CDR3HammingGraphComponentInfo hamming_graph_info(graph_component_map_,
                                                unique_cdr3s_map_,
                                                cdr3_to_old_index_map_,
                                                unique_cdr3s_,
                                                hg_component,
                                                component_id);
        std::shared_ptr<Base_CDR3_HG_CC_Processor> forest_calculator(
                new Kruskal_CDR3_HG_CC_Processor(clone_set_ptr_,
                                                 config_.algorithm_params,
                                                 clone_by_read_constructor_,
                                                 hamming_graph_info,
                                                 current_fake_clone_index_));
        auto tree = forest_calculator->Process();
        current_fake_clone_index_ = forest_calculator->GetCurrentFakeCloneIndex();
        reconstructed_ += forest_calculator->GetNumberOfReconstructedClones();
        return tree;
    }

    void VJClassProcessor::CreateUniqueCDR3Map() {
        const auto& clone_set = *clone_set_ptr_;
        for(auto it = decomposition_class_.begin(); it != decomposition_class_.end(); it++) {
            if(clone_set[*it].RegionIsEmpty(annotation_utils::StructuralRegion::CDR3))
                continue;
            auto cdr3 = core::dna5String_to_string(clone_set[*it].CDR3());
            if (unique_cdr3s_map_.find(cdr3) == unique_cdr3s_map_.end())
                unique_cdr3s_map_[cdr3] = std::vector<size_t>();
            unique_cdr3s_map_[cdr3].push_back(*it);
        }
        for(auto it = unique_cdr3s_map_.begin(); it != unique_cdr3s_map_.end(); it++)
            unique_cdr3s_.push_back(it->first);
        for(size_t i = 0; i < unique_cdr3s_.size(); ++i)
            cdr3_to_old_index_map_[unique_cdr3s_[i]] = i;
    }

    std::vector<SparseGraphPtr> VJClassProcessor::ComputeConnectedComponents() {
        CreateUniqueCDR3Map();
        std::string cdrs_fasta = WriteUniqueCDR3InFasta();
        std::string graph_fname = GetGraphFname();
        auto chain = BaseGeneClassProcessor::clone_set_ptr_->operator[](
                *decomposition_class_.cbegin()).ChainType().Chain();
        size_t tau = config_.algorithm_params.GetNumMismatchesByChainType(chain);
        return ComputeCDR3HammingGraphs(cdrs_fasta, graph_fname, tau);
    }




}