#include "v_class_processor.hpp"
#include <cdr3_hamming_graph_connected_components_processors/edmonds_cdr3_hg_cc_processor.hpp>


namespace antevolo {

    std::vector<SparseGraphPtr> VClassProcessor::ComputeConnectedComponents() {
        CreateUniqueCDR3Map();
        std::string cdrs_fasta = WriteUniqueCDR3InFasta();
        std::string graph_fname = GetGraphFname();

        auto chain = BaseGeneClassProcessor::clone_set_ptr_->operator[](
                *decomposition_class_.cbegin()).ChainType().Chain();
        size_t tau = config_.algorithm_params.GetNumMismatchesByChainType(chain);
        return ComputeCDR3HammingGraphs(cdrs_fasta, graph_fname, tau);
    }

    void VClassProcessor::CreateUniqueCDR3Map() {
        const auto &clone_set = *clone_set_ptr_;
        for (auto it = decomposition_class_.begin(); it != decomposition_class_.end(); it++) {
            if (clone_set[*it].RegionIsEmpty(annotation_utils::StructuralRegion::CDR3))
                continue;
            auto cdr3_jdifference_nucls = core::dna5String_to_string(
                    clone_set[*it].GetCDR3JDifferenceNucleotides(jdifference_positions_));
            if (unique_cdr3s_map_.find(cdr3_jdifference_nucls) == unique_cdr3s_map_.end())
                unique_cdr3s_map_[cdr3_jdifference_nucls] = std::vector<size_t>();
            unique_cdr3s_map_[cdr3_jdifference_nucls].push_back(*it);
        }
        for (auto it = unique_cdr3s_map_.begin(); it != unique_cdr3s_map_.end(); it++)
            unique_cdr3s_.push_back(it->first);
        for (size_t i = 0; i < unique_cdr3s_.size(); ++i)
            cdr3_to_old_index_map_[unique_cdr3s_[i]] = i;
    }

    void VClassProcessor::ChangeJgene(const germline_utils::ImmuneGene &v_gene,
                                      const germline_utils::ImmuneGene &j_gene) {
        auto &clone_set = *clone_set_ptr_;
        for (auto it = decomposition_class_.begin(); it != decomposition_class_.end(); it++) {
            if (clone_set[*it].RegionIsEmpty(annotation_utils::StructuralRegion::CDR3))
                continue;
            auto read = const_cast<core::Read &>(clone_set[*it].Read());
            clone_set[*it] = clone_by_read_constructor_.GetCloneByReadWithSpecificGenes(read, v_gene, j_gene);
        }
    }

    void VClassProcessor::ChangeJgeneToMax(CDR3HammingGraphComponentInfo hamming_graph_info) {
        auto &clone_set = *clone_set_ptr_;
        auto vertices = hamming_graph_info.GetAllClones();
        auto v_gene = clone_set[*vertices.begin()].VGene();
        // get most frequent J gene
        std::map<germline_utils::ImmuneGene, int> freq_map;
        germline_utils::ImmuneGene V;
        for (auto it = vertices.begin(); it != vertices.end(); it++) {
            auto it_find = freq_map.find(clone_set[*it].JGene());
            if (it_find == freq_map.end()) {
                freq_map.insert({clone_set[*it].JGene(), 0});
            }
            freq_map[clone_set[*it].JGene()]++;
            V = clone_set[*it].VGene();
        }
        germline_utils::ImmuneGene max_j_gene = std::max_element(
                std::begin(freq_map),
                std::end(freq_map),
                [](const std::pair<germline_utils::ImmuneGene, int> &p1,
                   const std::pair<germline_utils::ImmuneGene, int> &p2) {
                    return p1.second < p2.second;
                })->first;

        for (auto it = vertices.begin(); it != vertices.end(); it++) {
            if (clone_set[*it].RegionIsEmpty(annotation_utils::StructuralRegion::CDR3)) {
                continue;
            }

            auto &clone = clone_set[*it];
            std::string cdr3Jnucl = core::dna5String_to_string(
                    clone.GetCDR3JDifferenceNucleotides(jdifference_positions_));
            VERIFY(cdr3_to_old_index_map_.find(cdr3Jnucl) != cdr3_to_old_index_map_.end());
            auto read = const_cast<core::Read &>(clone_set[*it].Read());
            auto new_clone = clone_by_read_constructor_.GetCloneByReadWithSpecificGenes(read, v_gene, max_j_gene);
            clone_set[*it] = new_clone;
            std::string new_cdr3Jnucl = core::dna5String_to_string(
                    new_clone.GetCDR3JDifferenceNucleotides(jdifference_positions_));
            cdr3_to_old_index_map_.insert({new_cdr3Jnucl, cdr3_to_old_index_map_[cdr3Jnucl]});
            VERIFY(cdr3_to_old_index_map_.find(new_cdr3Jnucl) != cdr3_to_old_index_map_.end());
        }
    }

    EvolutionaryTree VClassProcessor::ProcessComponentWithEdmonds(SparseGraphPtr hg_component, size_t component_id,
                                                                  const ShmModelEdgeWeightCalculator &edge_weight_calculator) {

        CDR3HammingGraphComponentInfo hamming_graph_info(graph_component_map_,
                                                         unique_cdr3s_map_,
                                                         cdr3_to_old_index_map_,
                                                         unique_cdr3s_,
                                                         hg_component,
                                                         component_id);
        ChangeJgeneToMax(hamming_graph_info);
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


}