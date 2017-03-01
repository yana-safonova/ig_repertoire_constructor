#include "base_cdr3_hg_cc_processor.hpp"
#include <annotation_utils/shm_comparator.hpp>

namespace antevolo {

    Base_CDR3_HG_CC_Processor::Base_CDR3_HG_CC_Processor(
              CloneSetWithFakesPtr clone_set_ptr,
              const AntEvoloConfig::AlgorithmParams &config,
              const AnnotatedCloneByReadConstructor& clone_by_read_constructor,
              CDR3HammingGraphInfo& hamming_graph_info,
              size_t& current_fake_clone_index,
              size_t& reconstructed,
              size_t& rejected) :
            clone_set_ptr_(clone_set_ptr),
            config_(config),
            clone_by_read_constructor_(clone_by_read_constructor),
            hamming_graph_info_(hamming_graph_info),
            current_fake_clone_index_(current_fake_clone_index),
            reconstructed_(reconstructed),
            rejected_(rejected) { }

    EvolutionaryTree Base_CDR3_HG_CC_Processor::ConstructForest() {
        EvolutionaryTree tree(clone_set_ptr_);
        boost::unordered_set<size_t> vertices_nums(hamming_graph_info_.GetAllClones());


        size_t cdr3_length = clone_set_ptr_->operator[](*vertices_nums.cbegin()).CDR3Range().length();
        for (size_t clone_num : vertices_nums) {
            VERIFY(clone_set_ptr_->operator[](clone_num).CDR3Range().length() == cdr3_length);
        }


        std::map<size_t, size_t> rank;
        std::map<size_t, size_t> parent;
        boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges(rank, parent);
        for (size_t i : vertices_nums) {
            ds_on_undirected_edges.make_set(i);
        }
        AddUndirectedForest(ds_on_undirected_edges, vertices_nums);
        SetUndirectedComponentsParentEdges(ds_on_undirected_edges, vertices_nums);
        SetDirections(ds_on_undirected_edges, vertices_nums, tree);
        ReconstructMissingVertices(vertices_nums, tree);
        tree.AddAllEdges();
        return tree;
    }

    void Base_CDR3_HG_CC_Processor::AddUndirectedForest(boost::disjoint_sets<AP_map, AP_map> &ds_on_undirected_edges,
                                                        const boost::unordered_set<size_t> &vertices_nums) {
        const auto& clone_set = *clone_set_ptr_;
        for (auto src_num : vertices_nums) {
            size_t dst_num;
            auto it = getRelatedClonesIterator(hamming_graph_info_, clone_set[src_num]);
            while (it.HasNext()) {
                dst_num = it.Next();
                VERIFY(vertices_nums.find(dst_num) != vertices_nums.end());
                if (dst_num == src_num) {
                    continue;
                }
                if (ds_on_undirected_edges.find_set(src_num) !=
                    ds_on_undirected_edges.find_set(dst_num) &&
                    annotation_utils::SHMComparator::SHMsAreEqual(
                            clone_set[src_num].VSHMs(), clone_set[dst_num].VSHMs()) &&
                    annotation_utils::SHMComparator::SHMsAreEqual(
                            clone_set[src_num].JSHMs(), clone_set[dst_num].JSHMs())) {
                    AddUndirectedPair(src_num, dst_num);
                    ds_on_undirected_edges.union_set(src_num, dst_num);
                }
            }
        }
    }

    void Base_CDR3_HG_CC_Processor::AddUndirectedPair(size_t src_num, size_t dst_num) {
        if (undirected_graph_.find(src_num) == undirected_graph_.end()) {
            undirected_graph_[src_num] = std::set<size_t>();
        }
        if (undirected_graph_.find(dst_num) == undirected_graph_.end()) {
            undirected_graph_[dst_num] = std::set<size_t>();
        }
        if (undirected_graph_[src_num].find(dst_num) == undirected_graph_[src_num].end() &&
            undirected_graph_[dst_num].find(src_num) == undirected_graph_[dst_num].end()) {
            undirected_graph_[src_num].insert(dst_num);
            undirected_graph_[dst_num].insert(src_num);
        }
    }

}