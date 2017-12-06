#include "vj_class_processor.hpp"
#include <cdr3_hamming_graph_connected_components_processors/edmonds_cdr3_hg_cc_processor.hpp>


namespace antevolo {

    EvolutionaryTree VJClassProcessor::ProcessComponentWithKruskal(SparseGraphPtr hg_component, size_t component_id) {

        CDR3HammingGraphInfo hamming_graph_info(graph_component_map_,
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

    std::string VJClassProcessor::WriteUniqueCDR3InFasta(core::DecompositionClass decomposition_class) {
        std::string output_fname = GetFastaFname(decomposition_class);
        std::ofstream out(output_fname);
        for(auto it = unique_cdr3s_.begin(); it != unique_cdr3s_.end(); it++)
            out << ">" << *it << std::endl << *it << std::endl;
        return output_fname;
    }

    vector<SparseGraphPtr> VJClassProcessor::ComputeConnectedComponents(
            const core::DecompositionClass& decomposition_class) {
        CreateUniqueCDR3Map(decomposition_class);
        std::string cdrs_fasta = WriteUniqueCDR3InFasta(decomposition_class);
        std::string graph_fname = GetGraphFname(decomposition_class);
        return ComputeCDR3HammingGraphs(cdrs_fasta, graph_fname);
    }




}