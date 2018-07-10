#include "base_gene_class_processor.hpp"
#include <path_helper.hpp>
#include "../../graph_utils/graph_io.hpp"
#include "../../graph_utils/graph_splitter.hpp"
#include <convert.hpp>
#include <annotation_utils/shm_comparator.hpp>
#include <cdr3_hamming_graph_connected_components_processors/edmonds_cdr3_hg_cc_processor.hpp>

namespace antevolo {

    // return connected components of Hamming graph on CDR3s
    std::vector<SparseGraphPtr> BaseGeneClassProcessor::ComputeCDR3HammingGraphs(std::string cdr_fasta,
                                                                                 std::string graph_fname,
                                                                                 size_t tau) {
//        std::string run_graph_constructor = "./build/release/bin/ig_swgraph_construct";
        std::stringstream ss;
        ss << config_.cdr_labeler_config.input_params.run_hg_constructor << " -i " << cdr_fasta <<
           " -o " << graph_fname
           << " --tau " << tau
           << " --max-indels " << config_.algorithm_params.similar_cdr3s_params.num_indels
           << " -S " << " 0 " <<
           " -T " << " 0 " << " -k 10 > " << config_.output_params.trash_output;
        int err_code = system(ss.str().c_str());
        VERIFY_MSG(err_code == 0, "Graph constructor finished abnormally, error code: " << err_code);
        auto sparse_cdr_graph_ = GraphReader(graph_fname).CreateGraph();
        TRACE("Hamming graph contains " << sparse_cdr_graph_->N() << " edges and " << sparse_cdr_graph_->NZ()
                                        << " edges");

        auto connected_components = ConnectedComponentGraphSplitter(sparse_cdr_graph_).Split();
        graph_component_map_ = sparse_cdr_graph_->GetGraphComponentMap();
        return connected_components;
    }


    void BaseGeneClassProcessor::Clear() {
        unique_cdr3s_.clear();
        unique_cdr3s_map_.clear();
    }

    std::string BaseGeneClassProcessor::WriteUniqueCDR3InFasta() {
        std::string output_fname = GetFastaFname();
        std::ofstream out(output_fname);
        for (auto it = unique_cdr3s_.begin(); it != unique_cdr3s_.end(); it++)
            out << ">" << *it << std::endl << *it << std::endl;
        return output_fname;
    }

    std::string BaseGeneClassProcessor::GetFastaFname() {
        std::stringstream ss;
        size_t key = *decomposition_class_.begin();
        ss << "CDR3_sequences_key_" << key << ".fasta";
        return path::append_path(config_.output_params.cdr_graph_dir, ss.str());
    }

    std::string BaseGeneClassProcessor::GetGraphFname() {
        std::stringstream ss;
        size_t key = *decomposition_class_.begin();
        ss << "CDR3_sequences_key_" << key << ".graph";
        return path::append_path(config_.output_params.cdr_graph_dir, ss.str());
    }

    EvolutionaryTree BaseGeneClassProcessor::ProcessComponentWithEdmonds(
            SparseGraphPtr hg_component,
            size_t component_id,
            const ShmModelEdgeWeightCalculator &edge_weight_calculator) {

        CDR3HammingGraphComponentInfo hamming_graph_info(graph_component_map_,
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
        auto tree = forest_calculator->Process();
        current_fake_clone_index_ = forest_calculator->GetCurrentFakeCloneIndex();
        reconstructed_ += forest_calculator->GetNumberOfReconstructedClones();
        return tree;
    }

    EvolutionaryTree BaseGeneClassProcessor::ProcessComponent(
            SparseGraphPtr hg_component,
            size_t component_id,
            const ShmModelEdgeWeightCalculator &edge_weight_calculator) {
        return ProcessComponentWithEdmonds(hg_component, component_id, edge_weight_calculator);
    }
}
