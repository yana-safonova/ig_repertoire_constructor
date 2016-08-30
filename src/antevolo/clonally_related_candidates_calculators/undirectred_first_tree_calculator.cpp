#include <path_helper.hpp>
#include "undirectred_first_tree_calculator.hpp"
#include "../../graph_utils/graph_io.hpp"
#include "../../graph_utils/graph_splitter.hpp"
#include <convert.hpp>
#include <annotation_utils/shm_comparator.hpp>



namespace antevolo {
    void UndirectedFirstTreeCalculator::Clear() {
        unique_cdr3s_.clear();
        unique_cdr3s_map_.clear();
    }

    void UndirectedFirstTreeCalculator::CreateUniqueCDR3Map(
            core::DecompositionClass decomposition_class) {
        for(auto it = decomposition_class.begin(); it != decomposition_class.end(); it++) {
            if(clone_set_[*it].RegionIsEmpty(annotation_utils::StructuralRegion::CDR3))
                continue;
            auto cdr3 = core::dna5String_to_string(clone_set_[*it].CDR3());
            if(unique_cdr3s_map_.find(cdr3) == unique_cdr3s_map_.end())
                unique_cdr3s_map_[cdr3] = std::vector<size_t>();
            unique_cdr3s_map_[cdr3].push_back(*it);
        }
        for(auto it = unique_cdr3s_map_.begin(); it != unique_cdr3s_map_.end(); it++)
            unique_cdr3s_.push_back(it->first);
    }

    std::string UndirectedFirstTreeCalculator::GetFastaFname(core::DecompositionClass decomposition_class) {
        std::stringstream ss;
        size_t key = *decomposition_class.begin();
        ss << "CDR3_sequences_key_" << key << ".fasta";
        return path::append_path(output_params_.cdr_graph_dir, ss.str());
    }

    std::string UndirectedFirstTreeCalculator::GetGraphFname(core::DecompositionClass decomposition_class) {
        std::stringstream ss;
        size_t key = *decomposition_class.begin();
        ss << "CDR3_sequences_key_" << key << ".graph";
        return path::append_path(output_params_.cdr_graph_dir, ss.str());
    }

    std::string UndirectedFirstTreeCalculator::WriteUniqueCDR3InFasta(core::DecompositionClass decomposition_class) {
        std::string output_fname = GetFastaFname(decomposition_class);
        std::ofstream out(output_fname);
        for(auto it = unique_cdr3s_.begin(); it != unique_cdr3s_.end(); it++)
            out << ">" << *it << std::endl << *it << std::endl;
        return output_fname;
    }

    // return connected components of Hamming graph on CDR3s
    std::vector<SparseGraphPtr> UndirectedFirstTreeCalculator::ComputeCDR3HammingGraphs(std::string cdr_fasta,
                                                                                         std::string graph_fname) {
        std::string run_graph_constructor = "./build/release/bin/ig_swgraph_construct";
        std::stringstream ss;
        ss << run_graph_constructor << " -i " << cdr_fasta <<
                " -o " << graph_fname << " --tau " << num_mismatches_ << " -S " << " 0 " <<
                " -T " << " 0 " << " -k 10 > " << output_params_.trash_output;
        int err_code = system(ss.str().c_str());
        VERIFY_MSG(err_code == 0, "Graph constructor finished abnormally, error code: " << err_code);
        auto sparse_cdr_graph_ = GraphReader(graph_fname).CreateGraph();
        TRACE("Hamming graph contains " << sparse_cdr_graph_->N() << " edges and " << sparse_cdr_graph_->NZ() << " edges");
        //if(sparse_cdr_graph_->NZ() != 0)
        //    std::cout << *sparse_cdr_graph_ << std::endl;
        auto connected_components = ConnectedComponentGraphSplitter(sparse_cdr_graph_).Split();
        graph_component_ = sparse_cdr_graph_->GetGraphComponentMap();
        return connected_components;
    }


    EvolutionaryTree UndirectedFirstTreeCalculator::AddComponent(SparseGraphPtr hg_component,
                                                     size_t component_id) {

        auto forest_calculator = std::shared_ptr<BaseClusterToForestCalculator>(
                new KruskalClusterToForestCalculator(clone_set_,
                                                     config_,
                                                     graph_component_,
                                                     unique_cdr3s_map_,
                                                     unique_cdr3s_));
        auto tree = forest_calculator->ConstructForest(hg_component, component_id);
        return tree;
    }




}