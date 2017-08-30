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
        unique_cdr3sJnucl_map_.clear();
    }


//    void VJClassProcessor::CreateUniqueCDR3Map(
//            core::DecompositionClass decomposition_class) {
//        const auto &clone_set = *clone_set_ptr_;
//        for (auto it = decomposition_class.begin(); it != decomposition_class.end(); it++) {
//            if (clone_set[*it].RegionIsEmpty(annotation_utils::StructuralRegion::CDR3))
//                continue;
//            auto cdr3 = core::dna5String_to_string(clone_set[*it].CDR3());
//            if (unique_cdr3sJnucl_map_.find(cdr3) == unique_cdr3sJnucl_map_.end())
//                unique_cdr3sJnucl_map_[cdr3] = std::vector<size_t>();
//            unique_cdr3sJnucl_map_[cdr3].push_back(*it);
//        }
//        for (auto it = unique_cdr3sJnucl_map_.begin(); it != unique_cdr3sJnucl_map_.end(); it++)
//            unique_cdr3s_.push_back(it->first);
//        for (size_t i = 0; i < unique_cdr3s_.size(); ++i)
//            cdr3_to_old_index_map_[unique_cdr3s_[i]] = i;
//    }

    std::string VJClassProcessor::GetJNucleotides(const annotation_utils::AnnotatedClone &clone) {
        //indels! :(
        auto row_query = seqan::row(clone.JAlignment().Alignment(), 1);
        std::string JNucleotides;
        JNucleotides += core::dna5String_to_string(row_query[length(row_query) - 1 - 28]);
        JNucleotides += core::dna5String_to_string(row_query[length(row_query) - 1 - 27]);
        JNucleotides += core::dna5String_to_string(row_query[length(row_query) - 1 - 26]);
        JNucleotides += core::dna5String_to_string(row_query[length(row_query) - 1 - 25]);
        JNucleotides += core::dna5String_to_string(row_query[length(row_query) - 1 - 22]);
        JNucleotides += core::dna5String_to_string(row_query[length(row_query) - 1 - 19]);
        JNucleotides += core::dna5String_to_string(row_query[length(row_query) - 1 - 18]);
        JNucleotides += core::dna5String_to_string(row_query[length(row_query) - 1 - 17]);
        return JNucleotides;
    }

    void VJClassProcessor::CreateUniqueCDR3JNucleotidesMap(core::DecompositionClass decomposition_class) {
        const auto &clone_set = *clone_set_ptr_;
        for (auto it = decomposition_class.begin(); it != decomposition_class.end(); it++) {
            if (clone_set[*it].RegionIsEmpty(annotation_utils::StructuralRegion::CDR3))
                continue;
            auto cdr3 = core::dna5String_to_string(clone_set[*it].CDR3());
            auto JNucleotides = GetJNucleotides(clone_set[*it]);
            auto CDR3JNucleotides = std::make_pair(cdr3, JNucleotides);

            if (unique_cdr3sJnucl_map_.find(CDR3JNucleotides) == unique_cdr3sJnucl_map_.end())
                unique_cdr3sJnucl_map_[CDR3JNucleotides] = std::vector<size_t>();
            unique_cdr3sJnucl_map_[CDR3JNucleotides].push_back(*it);
        }
        for (auto it = unique_cdr3sJnucl_map_.begin(); it != unique_cdr3sJnucl_map_.end(); it++)
            unique_cdr3s_.push_back(it->first.first);
        for (size_t i = 0; i < unique_cdr3s_.size(); ++i)
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
        for (auto it = unique_cdr3s_.begin(); it != unique_cdr3s_.end(); it++)
            out << ">" << *it << std::endl << *it << std::endl;
        return output_fname;
    }


    std::string WriteUniqueCDR3JNucleotidesInFasta(core::DecompositionClass decomposition_class);


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
        TRACE("Hamming graph contains " << sparse_cdr_graph_->N() << " edges and " << sparse_cdr_graph_->NZ()
                                        << " edges");

        auto connected_components = ConnectedComponentGraphSplitter(sparse_cdr_graph_).Split();
        graph_component_map_ = sparse_cdr_graph_->GetGraphComponentMap();
        return connected_components;
    }


    double VJClassProcessor::CDR3JNucleotidesDistance(const std::pair<string, string> &cdr3_jnucl1,
                                                      const std::pair<string, string> &cdr3_jnucl2) {
        int jnucl_mismatches = 0;
        for (int i = 0; i < cdr3_jnucl1.second.size(); ++i) {
            if ((cdr3_jnucl1.second[i] != cdr3_jnucl2.second[i]) ||
                (cdr3_jnucl1.second[i] == 'N' && cdr3_jnucl2.second[i] == 'N'))
                ++jnucl_mismatches;
        }
        return jnucl_mismatches;
    }

    void VJClassProcessor::dfs(int v, bool *used, std::vector<int> &components_helper,
                                   std::vector<std::pair<string, string>> const &unique_cdr3Jnucl, double tau) {
        used[v] = true;
        components_helper.push_back(v);
        for (size_t i = 0; i < unique_cdr3Jnucl.size(); ++i) {
            if (!used[i] && CDR3JNucleotidesDistance(unique_cdr3Jnucl[v], unique_cdr3Jnucl[i]) < tau)
                dfs(i, used, components_helper, unique_cdr3Jnucl, tau);
        }
    }


// return connected components of Hamming graph on CDR3s + JNucleotides
    void VJClassProcessor::ComputeCDR3JNucleotidesHammingGraphs(double tau, std::string graph_fname) {
        std::vector<std::pair<string, string>> unique_cdr3Jnucl;
        for (auto const &it:unique_cdr3sJnucl_map_)
            unique_cdr3Jnucl.push_back(it.first);
        std::map<std::vector<std::pair<string, string>>, std::vector<int>> connected_components;
        std::vector<int> components_helper;

        bool *used = new bool[unique_cdr3Jnucl.size()];
        for (int i = 0; i < unique_cdr3Jnucl.size(); ++i)
            if (!used[i]) {
                components_helper.clear();
                dfs(i, used, components_helper, unique_cdr3Jnucl, tau);

//                cout << "Component:";
//                for (size_t j = 0; j < comp.size(); ++j)
//                    cout << ' ' << comp[j];
//                cout << endl;
            }
        delete[] used;





//        m.lock();
//        for (auto it = unique_cdr3sJnucl_map_.begin(); it!=unique_cdr3sJnucl_map_.end(); ++it)
//            if (it->second.size() > 1)
//                std::cout << it->first.first << " " << it->first.second << " => " << it->second.size() << '\n';
//        m.unlock();
////        std::string run_graph_constructor = "./build/release/bin/ig_swgraph_construct";
//        std::stringstream ss;
//        ss << config_.cdr_labeler_config.input_params.run_hg_constructor << " -i " << cdr_fasta <<
//           " -o " << graph_fname << " --tau " << num_mismatches_ << " -S " << " 0 " <<
//           " -T " << " 0 " << " -k 10 > " << config_.output_params.trash_output;
//        int err_code = system(ss.str().c_str());
//        VERIFY_MSG(err_code == 0, "Graph constructor finished abnormally, error code: " << err_code);
//        auto sparse_cdr_graph_ = GraphReader(graph_fname).CreateGraph();
//        TRACE("Hamming graph contains " << sparse_cdr_graph_->N() << " edges and " << sparse_cdr_graph_->NZ()
//                                        << " edges");
//
//        auto connected_components = ConnectedComponentGraphSplitter(sparse_cdr_graph_).Split();
//        graph_component_map_ = sparse_cdr_graph_->GetGraphComponentMap();
//        return connected_components;
    }


    void VJClassProcessor::ProcessComponentWithKruskal(SparseGraphPtr hg_component, size_t component_id) {

//        CDR3HammingGraphInfo hamming_graph_info(graph_component_map_,
//                                                unique_cdr3sJnucl_map_,
//                                                cdr3_to_old_index_map_,
//                                                unique_cdr3s_,
//                                                hg_component,
//                                                component_id);
//        std::shared_ptr<Base_CDR3_HG_CC_Processor> forest_calculator(
//                new Kruskal_CDR3_HG_CC_Processor(clone_set_ptr_,
//                                                 config_.algorithm_params,
//                                                 clone_by_read_constructor_,
//                                                 hamming_graph_info,
//                                                 current_fake_clone_index_));
//        auto tree = forest_calculator->ConstructForest();
//        current_fake_clone_index_ = forest_calculator->GetCurrentFakeCloneIndex();
//        reconstructed_ += forest_calculator->GetNumberOfReconstructedClones();
//        return tree;
    }

    void VJClassProcessor::ProcessComponentWithEdmonds(
            //EvolutionaryTree VJClassProcessor::ProcessComponentWithEdmonds(
            SparseGraphPtr hg_component,
            size_t component_id,
            const ShmModelEdgeWeightCalculator &edge_weight_calculator) {

//        CDR3HammingGraphInfo hamming_graph_info(graph_component_map_,
//                                                unique_cdr3sJnucl_map_,
//                                                cdr3_to_old_index_map_,
//                                                unique_cdr3s_,
//                                                hg_component,
//                                                component_id);
//        std::shared_ptr<Base_CDR3_HG_CC_Processor> forest_calculator(
//                new Edmonds_CDR3_HG_CC_Processor(clone_set_ptr_,
//                                                 config_.algorithm_params,
//                                                 clone_by_read_constructor_,
//                                                 hamming_graph_info,
//                                                 current_fake_clone_index_,
//                                                 edge_weight_calculator));
//        auto tree = forest_calculator->ConstructForest();
//        current_fake_clone_index_ = forest_calculator->GetCurrentFakeCloneIndex();
//        reconstructed_ += forest_calculator->GetNumberOfReconstructedClones();
//        return tree;
    }


    void VJClassProcessor::HG_components(SparseGraphPtr hg_component, size_t component_id,
                                         const ShmModelEdgeWeightCalculator &edge_weight_calculator) {
//        CDR3HammingGraphInfo hamming_graph_info(graph_component_map_,
//                                                unique_cdr3sJnucl_map_,
//                                                cdr3_to_old_index_map_,
//                                                unique_cdr3s_,
//                                                hg_component,
//                                                component_id);
//
//        auto clones = hamming_graph_info.GetAllClones();
//        const auto &clone_set = *clone_set_ptr_;
//
//        std::string output_fname = path::append_path(config_.output_params.output_dir, "HG_stats.txt");
//        std::ofstream out(output_fname, std::ios::app);
//
//        out << "Component size: " << clones.size() << "\n";
//        std::map<std::string, int> q;
//        size_t cdr3_len;
//        for (auto const &it: clones) {
//            std::string v_name = std::string(seqan::toCString(clone_set[it].VGene().name()));
//            std::string j_name = std::string(seqan::toCString(clone_set[it].JGene().name()));
//            q[v_name + '_' + j_name] += 1;
//            cdr3_len = seqan::length(clone_set[it].CDR3());
//        }
//        for (auto it = q.begin(); it != q.end(); ++it) {
//            out << it->first << " " << it->second << "\n";
//        }
//        out << "CDR3: " << cdr3_len << "\n\n";
//        out.close();
    }


}