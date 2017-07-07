#include <path_helper.hpp>
#include "vj_class_processor.hpp"
#include "../../graph_utils/graph_io.hpp"
#include "../../graph_utils/graph_splitter.hpp"
#include <convert.hpp>
#include <annotation_utils/shm_comparator.hpp>
#include <cdr3_hamming_graph_connected_components_processors/edmonds_cdr3_hg_cc_processor.hpp>


namespace antevolo {
    void VJClassProcessor::Clear() {
        unique_cdr3s_.clear();
        unique_cdr3s_map_.clear();
    }

    void VJClassProcessor::CreateUniqueCDR3Map(
            core::DecompositionClass decomposition_class) {
        const auto& clone_set = *clone_set_ptr_;
        for(auto it = decomposition_class.begin(); it != decomposition_class.end(); it++) {
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

    std::string VJClassProcessor::GetFastaFname(core::DecompositionClass decomposition_class) {
        std::stringstream ss;
        size_t key = *decomposition_class.begin();
        ss << "CDR3_sequences_key_" << key << ".fasta";
        return path::append_path(config_.output_params.cdr_graph_dir, ss.str());
    }

    std::string VJClassProcessor::GetGraphFname(core::DecompositionClass decomposition_class) {
        std::stringstream ss;
        size_t key = *decomposition_class.begin();
        ss << "CDR3_sequences_key_" << key << ".graph";
        return path::append_path(config_.output_params.cdr_graph_dir, ss.str());
    }

    std::string VJClassProcessor::WriteUniqueCDR3InFasta(core::DecompositionClass decomposition_class) {
        std::string output_fname = GetFastaFname(decomposition_class);
        std::ofstream out(output_fname);
        for(auto it = unique_cdr3s_.begin(); it != unique_cdr3s_.end(); it++)
            out << ">" << *it << std::endl << *it << std::endl;
        return output_fname;
    }

    // return connected components of Hamming graph on CDR3s
    std::vector<SparseGraphPtr> VJClassProcessor::ComputeCDR3HammingGraphs(std::string cdr_fasta,
                                                                           std::string graph_fname) {
//        std::string run_graph_constructor = "./build/release/bin/ig_swgraph_construct";
        std::stringstream ss;
        ss << config_.cdr_labeler_config.input_params.run_hg_constructor << " -i " << cdr_fasta <<
                " -o " << graph_fname << " --tau " << num_mismatches_ << " -S " << " 0 " <<
                " -T " << " 0 " << " -k 10 > " << config_.output_params.trash_output;
        int err_code = system(ss.str().c_str());
        VERIFY_MSG(err_code == 0, "Graph constructor finished abnormally, error code: " << err_code);
        auto sparse_cdr_graph_ = GraphReader(graph_fname).CreateGraph();
        TRACE("Hamming graph contains " << sparse_cdr_graph_->N() << " edges and " << sparse_cdr_graph_->NZ() << " edges");

        auto connected_components = ConnectedComponentGraphSplitter(sparse_cdr_graph_).Split();
        graph_component_map_ = sparse_cdr_graph_->GetGraphComponentMap();
        return connected_components;
    }


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
        auto tree = forest_calculator->ConstructForest();
        current_fake_clone_index_ = forest_calculator->GetCurrentFakeCloneIndex();
        reconstructed_ += forest_calculator->GetNumberOfReconstructedClones();
        return tree;
    }
    EvolutionaryTree VJClassProcessor::ProcessComponentWithEdmonds(
            SparseGraphPtr hg_component,
            size_t component_id,
            const ShmModelEdgeWeightCalculator& edge_weight_calculator) {

        CDR3HammingGraphInfo hamming_graph_info(graph_component_map_,
                                                unique_cdr3s_map_,
                                                cdr3_to_old_index_map_,
                                                unique_cdr3s_,
                                                hg_component,
                                                component_id);
        std::shared_ptr<Base_CDR3_HG_CC_Processor> forest_calculator(
                    new Edmonds_CDR3_HG_CC_Processor(clone_set_ptr_,
                                                     config_.algorithm_params,
                                                     clone_by_read_constructor_,
                                                     hamming_graph_info,
                                                     current_fake_clone_index_,
                                                     edge_weight_calculator));
        auto tree = forest_calculator->ConstructForest();
        current_fake_clone_index_ = forest_calculator->GetCurrentFakeCloneIndex();
        reconstructed_ += forest_calculator->GetNumberOfReconstructedClones();
        return tree;
    }




}