#include "base_cdr3_hg_cc_processor.hpp"
#include <annotation_utils/shm_comparator.hpp>

namespace antevolo {

    Base_CDR3_HG_CC_Processor::Base_CDR3_HG_CC_Processor(
                                      CloneSetWithFakes &clone_set,
                                      const AntEvoloConfig::AlgorithmParams &config,
                                      const AnnotatedCloneByReadConstructor& clone_by_read_constructor,
                                      CDR3HammingGraphInfo& hamming_graph_info) :
            clone_set_(clone_set),
            config_(config),
            clone_by_read_constructor_(clone_by_read_constructor),
            hamming_graph_info_(hamming_graph_info) { }

    EvolutionaryTree Base_CDR3_HG_CC_Processor::ConstructForest() {
        EvolutionaryTree tree(clone_set_);
        boost::unordered_set<size_t> vertices_nums(hamming_graph_info_.GetAllClones());
        std::map<size_t, size_t> rank;
        std::map<size_t, size_t> parent;
        boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges(rank, parent);
        for (size_t i : vertices_nums) {
            ds_on_undirected_edges.make_set(i);
        }
        AddUndirectedForest(ds_on_undirected_edges, vertices_nums);
        SetUndirectedComponentsParentEdges(ds_on_undirected_edges, vertices_nums);
        SetDirections(ds_on_undirected_edges, vertices_nums, tree);
//        ReconstructMissingVertices(vertices_nums, tree, hg_component, component_id);
        tree.AddAllEdges();
        return tree;
    }

    void Base_CDR3_HG_CC_Processor::AddUndirectedForest(boost::disjoint_sets<AP_map, AP_map>& ds_on_undirected_edges,
                                                        boost::unordered_set<size_t> vertices_nums) {
        for (auto src_num : vertices_nums) {
            size_t dst_num;
            auto it = getRelatedClonesIterator(hamming_graph_info_, clone_set_[src_num]);
            while (it.HasNext()) {
                dst_num = it.Next();
                VERIFY(vertices_nums.find(dst_num) != vertices_nums.end());
                if (dst_num == src_num) {
                    continue;
                }
                if (ds_on_undirected_edges.find_set(src_num) !=
                    ds_on_undirected_edges.find_set(dst_num) &&
                    annotation_utils::SHMComparator::SHMsAreEqual(
                            clone_set_[src_num].VSHMs(), clone_set_[dst_num].VSHMs()) &&
                    annotation_utils::SHMComparator::SHMsAreEqual(
                            clone_set_[src_num].JSHMs(), clone_set_[dst_num].JSHMs())) {
                    AddUndirectedPair(src_num, dst_num);
                    ds_on_undirected_edges.union_set(src_num, dst_num);
                }
            }
        }
        /*
        //adding undirected edges first
        for (size_t i = 0; i < hg_component->N(); i++) {
            size_t old_index = graph_component_map_.GetOldVertexByNewVertex(component_id, i);
            //auto clones_sharing_cdr3 = unique_cdr3s_map_[unique_cdr3s_[old_index]];
            auto clones_sharing_cdr3 = unique_cdr3s_map_.find(unique_cdr3s_[old_index])->second;
            for (size_t it1 = 0; it1 < clones_sharing_cdr3.size(); it1++) {
                parent_edge_handled_[clones_sharing_cdr3[it1]] = false;
                for (size_t it2 = it1 + 1; it2 < clones_sharing_cdr3.size(); it2++) {
                    size_t src_num = clones_sharing_cdr3[it1];
                    size_t dst_num = clones_sharing_cdr3[it2];
                    // if clones are not in the same connected component and
                    // the edge is undirected
                    if (ds_on_undirected_edges.find_set(src_num) !=
                        ds_on_undirected_edges.find_set(dst_num) &&
                        annotation_utils::SHMComparator::SHMsAreEqual(
                                clone_set_[src_num].VSHMs(), clone_set_[dst_num].VSHMs()) &&
                        annotation_utils::SHMComparator::SHMsAreEqual(
                                clone_set_[src_num].JSHMs(), clone_set_[dst_num].JSHMs())) {
                        AddUndirectedPair(src_num, dst_num);
                        ds_on_undirected_edges.union_set(src_num, dst_num);
                    }
                }
            }
        }
        for(size_t i = 0; i < hg_component->N(); i++) {
            size_t old_index1 = graph_component_map_.GetOldVertexByNewVertex(component_id, i);
            for (auto it = hg_component->VertexEdges(i).begin(); it != hg_component->VertexEdges(i).end(); it++) {
                size_t old_index2 = graph_component_map_.GetOldVertexByNewVertex(component_id, *it);
                //auto indices_1 = unique_cdr3s_map_[unique_cdr3s_[old_index1]];
                auto indices_1 = unique_cdr3s_map_.find(unique_cdr3s_[old_index1])->second;
                //auto indices_2 = unique_cdr3s_map_[unique_cdr3s_[old_index2]];
                auto indices_2 = unique_cdr3s_map_.find(unique_cdr3s_[old_index2])->second;
                for (auto it1 = indices_1.begin(); it1 != indices_1.end(); it1++)
                    for (auto it2 = indices_2.begin(); it2 != indices_2.end(); it2++) {
                        if (ds_on_undirected_edges.find_set(*it1) !=
                            ds_on_undirected_edges.find_set(*it2) &&
                            annotation_utils::SHMComparator::SHMsAreEqual(
                                    clone_set_[*it1].VSHMs(),
                                    clone_set_[*it2].VSHMs()) &&
                            annotation_utils::SHMComparator::SHMsAreEqual(
                                    clone_set_[*it1].JSHMs(),
                                    clone_set_[*it2].JSHMs())) {
                            AddUndirectedPair(*it1, *it2);
                            ds_on_undirected_edges.union_set(*it1, *it2);
                        }
                    }
            }
        }
        */
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