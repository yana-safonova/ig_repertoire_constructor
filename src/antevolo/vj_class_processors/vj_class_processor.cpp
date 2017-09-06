#include <path_helper.hpp>
#include "vj_class_processor.hpp"
#include "../../graph_utils/graph_io.hpp"
#include "../../graph_utils/graph_splitter.hpp"
#include <convert.hpp>
#include <annotation_utils/shm_comparator.hpp>
#include <cdr3_hamming_graph_connected_components_processors/edmonds_cdr3_hg_cc_processor.hpp>
#include <seqan/align.h>
#include <mutex>

namespace antevolo {
    void VJClassProcessor::Clear() {
        unique_cdr3Jnucleotides_.clear();
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
//            unique_cdr3Jnucleotides_.push_back(it->first);
//        for (size_t i = 0; i < unique_cdr3Jnucleotides_.size(); ++i)
//            cdr3Jnucl_to_old_index_map_[unique_cdr3Jnucleotides_[i]] = i;
//    }

    std::mutex m;

    seqan::Dna5String VJClassProcessor::GetJNucleotides(const annotation_utils::AnnotatedClone &clone) {
        //indels! :(
        auto row_query = seqan::row(clone.JAlignment().Alignment(), 1);
//        std::cout << row_query << "\n";
//        m.unlock();
        std::string JNucleotides(8, 'N');
        JNucleotides[0] = core::dna5String_to_string(row_query[length(row_query) - 1 - 28])[0];
        JNucleotides[1] = core::dna5String_to_string(row_query[length(row_query) - 1 - 27])[0];
        JNucleotides[2] = core::dna5String_to_string(row_query[length(row_query) - 1 - 26])[0];
        JNucleotides[3] = core::dna5String_to_string(row_query[length(row_query) - 1 - 25])[0];
        JNucleotides[4] = core::dna5String_to_string(row_query[length(row_query) - 1 - 22])[0];
        JNucleotides[5] = core::dna5String_to_string(row_query[length(row_query) - 1 - 19])[0];
        JNucleotides[6] = core::dna5String_to_string(row_query[length(row_query) - 1 - 18])[0];
        JNucleotides[7] = core::dna5String_to_string(row_query[length(row_query) - 1 - 17])[0];
        return seqan::Dna5String(JNucleotides);
    }

    void VJClassProcessor::CreateUniqueCDR3JNucleotidesMap(core::DecompositionClass decomposition_class) {
        const auto &clone_set = *clone_set_ptr_;
        for (auto it = decomposition_class.begin(); it != decomposition_class.end(); it++) {
            if (clone_set[*it].RegionIsEmpty(annotation_utils::StructuralRegion::CDR3))
                continue;
            auto cdr3 = clone_set[*it].CDR3();
            auto JNucleotides = GetJNucleotides(clone_set[*it]);
            auto CDR3JNucleotides = std::make_pair(cdr3, JNucleotides);

            if (unique_cdr3sJnucl_map_.find(CDR3JNucleotides) == unique_cdr3sJnucl_map_.end())
                unique_cdr3sJnucl_map_[CDR3JNucleotides] = std::vector<size_t>();
            unique_cdr3sJnucl_map_[CDR3JNucleotides].push_back(*it);
        }

        for (auto it = unique_cdr3sJnucl_map_.begin(); it != unique_cdr3sJnucl_map_.end(); it++)
            unique_cdr3Jnucleotides_.push_back(it->first);
        for (size_t i = 0; i < unique_cdr3Jnucleotides_.size(); ++i)
            cdr3Jnucl_to_old_index_map_[unique_cdr3Jnucleotides_[i]] = i;
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
        for (auto it = unique_cdr3Jnucleotides_.begin(); it != unique_cdr3Jnucleotides_.end(); it++)
            out << ">" << *it << std::endl << *it << std::endl;
        return output_fname;
    }


    std::string VJClassProcessor::WriteUniqueCDR3JNucleotidesInFasta(core::DecompositionClass decomposition_class)
    {
        std::string output_fname = GetFastaFname(decomposition_class);
        std::ofstream out(output_fname);
        for (auto it = unique_cdr3Jnucleotides_.begin(); it != unique_cdr3Jnucleotides_.end(); it++)
            out << ">" << it->first << it->second << "\n" << it->first << it->second << "\n";
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
        TRACE("Hamming graph contains " << sparse_cdr_graph_->N() << " edges and " << sparse_cdr_graph_->NZ()
                                        << " edges");

        auto connected_components = ConnectedComponentGraphSplitter(sparse_cdr_graph_).Split();
        graph_component_map_ = sparse_cdr_graph_->GetGraphComponentMap();
        return connected_components;
    }


    double
    VJClassProcessor::CDR3JNucleotidesDistance(const std::pair<seqan::Dna5String, seqan::Dna5String> &cdr3_jnucl1,
                                               const std::pair<seqan::Dna5String, seqan::Dna5String> &cdr3_jnucl2) {
        //m.lock();

        //CDR3 levenshtein distance
        seqan::Score<int> score(0, -1, -1, -1); //  match = 0, mismatch = -1, gapextend = -1, gapopen = -1
        seqan::Align<seqan::Dna5String, seqan::ArrayGaps> align;
        seqan::resize(seqan::rows(align), 2);
        seqan::assignSource(seqan::row(align, 0), cdr3_jnucl1.first);
        seqan::assignSource(seqan::row(align, 1), cdr3_jnucl2.first);
        int cdr3_dist = -seqan::globalAlignment(align, score);

        // J nucleotides distance
        int jnucl_mismatches = 0;
        for (int i = 0; i < length(cdr3_jnucl1.second); ++i) {
            if ((cdr3_jnucl1.second[i] != cdr3_jnucl2.second[i]) ||
                (cdr3_jnucl1.second[i] == 'N' && cdr3_jnucl2.second[i] == 'N'))
                ++jnucl_mismatches;
        }

//        if (cdr3_dist < 10)
//            std::cout << cdr3_jnucl1.first << "\n" << cdr3_jnucl2.first << "\n" << cdr3_dist << "\n";
//        //std::cout << cdr3_jnucl1.second << "\n" << cdr3_jnucl2.second << "\n" << jnucl_mismatches << "\n\n\n";
        //m.unlock();

        return jnucl_mismatches + cdr3_dist; // maybe some f(jnucl_mismatches, cdr3_dist)
    }

    void VJClassProcessor::dfs(int v, bool *used, std::vector<int> &components_helper,
                               std::vector<std::pair<seqan::Dna5String, seqan::Dna5String>> const &unique_cdr3Jnucl,
                               double tau) {
        used[v] = true;
        components_helper.push_back(v);
        for (size_t i = 0; i < unique_cdr3Jnucl.size(); ++i) {
            if (!used[i] && CDR3JNucleotidesDistance(unique_cdr3Jnucl[v], unique_cdr3Jnucl[i]) < tau)
                dfs(i, used, components_helper, unique_cdr3Jnucl, tau);
        }
    }


// return connected components of Hamming graph on CDR3s + JNucleotides
    map<vector<std::pair<seqan::Dna5String, seqan::Dna5String>>, vector<int>>
    VJClassProcessor::ComputeCDR3JNucleotidesHammingGraphs(
            double tau, std::string graph_fname) {
        std::vector<std::pair<seqan::Dna5String, seqan::Dna5String>> unique_cdr3Jnucl;
        for (auto const &it:unique_cdr3sJnucl_map_)
            unique_cdr3Jnucl.push_back(it.first);
        std::map<std::vector<std::pair<seqan::Dna5String, seqan::Dna5String>>, std::vector<int>> connected_components;
        std::vector<int> components_helper_ind;

        bool *used = new bool[unique_cdr3Jnucl.size()];
        for (int i = 0; i < unique_cdr3Jnucl.size(); ++i)
            if (!used[i]) {
                components_helper_ind.clear();
                dfs(i, used, components_helper_ind, unique_cdr3Jnucl, tau);

                std::vector<std::pair<seqan::Dna5String, seqan::Dna5String>> components_key;
                std::vector<int> components_value;

                for (auto const &it:components_helper_ind) {
                    components_key.push_back(unique_cdr3Jnucl[it]);
                    components_value.insert(components_value.end(),
                                            unique_cdr3sJnucl_map_[unique_cdr3Jnucl[it]].begin(),
                                            unique_cdr3sJnucl_map_[unique_cdr3Jnucl[it]].end());
                }

                std::cout << components_key.size() << " " << components_value.size() << "\n";
                connected_components[components_key] = components_value;
            }
        delete[] used;
        return connected_components;
    }


    void VJClassProcessor::ProcessComponentWithKruskal(SparseGraphPtr hg_component, size_t component_id) {

//        CDR3HammingGraphInfo hamming_graph_info(graph_component_map_,
//                                                unique_cdr3sJnucl_map_,
//                                                cdr3Jnucl_to_old_index_map_,
//                                                unique_cdr3Jnucleotides_,
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

    //void VJClassProcessor::ProcessComponentWithEdmonds(
    EvolutionaryTree VJClassProcessor::ProcessComponentWithEdmonds(
            SparseGraphPtr hg_component,
            size_t component_id,
            const ShmModelEdgeWeightCalculator &edge_weight_calculator) {

        CDR3HammingGraphInfo hamming_graph_info(graph_component_map_,
                                                unique_cdr3sJnucl_map_,
                                                cdr3Jnucl_to_old_index_map_,
                                                unique_cdr3Jnucleotides_,
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


    void VJClassProcessor::HG_components(SparseGraphPtr hg_component, size_t component_id,
                                         const ShmModelEdgeWeightCalculator &edge_weight_calculator) {
//        CDR3HammingGraphInfo hamming_graph_info(graph_component_map_,
//                                                unique_cdr3sJnucl_map_,
//                                                cdr3Jnucl_to_old_index_map_,
//                                                unique_cdr3Jnucleotides_,
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