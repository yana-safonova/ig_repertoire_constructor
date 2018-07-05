#include "base_cdr3_hg_cc_processor.hpp"
#include <annotation_utils/shm_comparator.hpp>
#include "parent_read_reconstructor.hpp"

namespace antevolo {
    Base_CDR3_HG_CC_Processor::Base_CDR3_HG_CC_Processor(
            CloneSetWithFakesPtr clone_set_ptr,
            const AntEvoloConfig::AlgorithmParams &config,
            const AnnotatedCloneByReadConstructor &clone_by_read_constructor,
            CDR3HammingGraphComponentInfo &hamming_graph_info,
            size_t current_fake_clone_index) :
            clone_set_ptr_(clone_set_ptr),
            config_(config),
            clone_by_read_constructor_(clone_by_read_constructor),
            hamming_graph_info_(hamming_graph_info),
            current_fake_clone_index_(current_fake_clone_index),
            reconstructed_(0) {}


    void Base_CDR3_HG_CC_Processor::AddUndirectedForest(boost::disjoint_sets<AP_map, AP_map> &ds_on_undirected_edges,
                                                        const boost::unordered_set<size_t> &vertices_nums) {
        const auto &clone_set = *clone_set_ptr_;
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

    void Base_CDR3_HG_CC_Processor::ReconstructMissingVertices(boost::unordered_set<size_t> &vertices_nums,
                                                               EvolutionaryTree &tree) {
        const auto &clone_set = *clone_set_ptr_;
        auto edge_constructor = GetEdgeConstructor();
        boost::unordered_map<size_t, EvolutionaryEdgePtr> roots_nearest_neighbours;
        std::vector<size_t> roots;
        boost::unordered_map<size_t, size_t> iterator_index_map;
        boost::unordered_set<size_t> rejected_roots;
        for (auto v : vertices_nums) {
            if (!tree.HasParentEdge(v)) {
                roots.push_back(v);
                iterator_index_map[v] = v;
            }
        }

        for (size_t root_num : roots) {
//            auto it = getRelatedClonesIterator(hamming_graph_info_, clone_set[root_num]);
//            while (it.HasNext()) {
//                size_t dst_num = it.Next();
//                VERIFY(vertices_nums.find(dst_num) != vertices_nums.end());
            for (size_t dst_num : vertices_nums) {
                if (dst_num == root_num) {
                    continue;
                }
                auto edge = edge_constructor->ConstructEdge(
                        clone_set[root_num],
                        clone_set[dst_num],
                        root_num,
                        dst_num);
                if (edge->IsIntersected()) {
                    if (roots_nearest_neighbours.find(root_num) == roots_nearest_neighbours.end() ||
                        edge->Length() < roots_nearest_neighbours[root_num]->Length()) {
                        roots_nearest_neighbours[root_num] = edge;
                    }
                }
            }
        }
        while (true) {
            size_t best_root_index = 0;
            size_t best_root_edge_length = EVO_EDGE_MAX_LENGTH;
            for (size_t i = 0; i < roots.size(); ++i) {
                if (!tree.HasParentEdge(roots[i]) &&
                    rejected_roots.find(roots[i]) == rejected_roots.end() &&
                    roots_nearest_neighbours.find(roots[i]) != roots_nearest_neighbours.end() &&
                    roots_nearest_neighbours[roots[i]]->Length() < best_root_edge_length) {
                    best_root_index = i;
                    best_root_edge_length = roots_nearest_neighbours[roots[i]]->Length();
                }
            }
            if (best_root_edge_length == EVO_EDGE_MAX_LENGTH) {
                break;
            }
            size_t root_num = roots[best_root_index];
            auto edge = roots_nearest_neighbours[root_num];

            if (rejected_roots.find(root_num) != rejected_roots.end()) {
                continue;
            }
            bool reconstructed = ReconstructAncestralLineageSimple(edge,
                                                                   tree,
                                                                   vertices_nums,
                                                                   edge_constructor,
                                                                   roots,
                                                                   iterator_index_map);
            if (!reconstructed) {
                rejected_roots.insert(root_num);
                continue;
            }
            root_num = clone_set.size() - 1;
            VERIFY(vertices_nums.find(root_num) != vertices_nums.end());
            for (size_t dst_num : vertices_nums) {
                HandleRootNeighbour(root_num,
                                    dst_num,
                                    vertices_nums,
                                    tree,
                                    roots_nearest_neighbours,
                                    edge_constructor);
            }
        }
    }

    void Base_CDR3_HG_CC_Processor::Refine(boost::unordered_set<size_t> &vertices_nums,
                                           EvolutionaryTree &tree) {
        auto &clone_set = *clone_set_ptr_;
        auto edge_constructor = GetEdgeConstructor();
        boost::unordered_map<size_t, EvolutionaryEdgePtr> best_intersected_edges;
        boost::unordered_map<size_t, EvolutionaryEdgePtr> best_reverse_edges;
//        INFO("start refinement");
        for (size_t clone_num : vertices_nums) {
//            auto it = getRelatedClonesIterator(hamming_graph_info_, clone_set[clone_num]);
//            while (it.HasNext()) {
//                size_t src_num = it.Next();
            for (size_t src_num : vertices_nums) {
                VERIFY(vertices_nums.find(src_num) != vertices_nums.end());
                if (src_num == clone_num) {
                    continue;
                }
                auto edge = edge_constructor->ConstructEdge(
                        clone_set[src_num],
                        clone_set[clone_num],
                        src_num,
                        clone_num);

                auto chain = clone_set[clone_num].ChainType().Chain();
                if (edge->CDR3Distance() > config_.GetNumMismatchesByChainType(chain) +
                                           config_.similar_cdr3s_params.num_indels) {
                    continue;
                }
                if (edge->IsIntersected()) {
                    if ((best_intersected_edges.find(clone_num) == best_intersected_edges.end() ||
                         edge->Length() < best_intersected_edges[clone_num]->Length()) &&
                        CheckClonesConsistencyForReconstruction(*edge->SrcClone(), *edge->DstClone())) {
                        best_intersected_edges[clone_num] = edge;
                    }
                }
                if (edge->IsReverseDirected()) {
                    if (best_reverse_edges.find(clone_num) == best_reverse_edges.end() ||
                        edge->Length() < best_reverse_edges[clone_num]->Length()) {
                        best_reverse_edges[clone_num] = edge;
                    }
                }
            }
        }
//        INFO("calculated best edges");
        std::vector<size_t> vertices_list;
        for (size_t clone_num : vertices_nums) {
            if (clone_set.IsFake(clone_num)) {
                continue;
            }
            vertices_list.push_back(clone_num);
        }
        for (size_t clone_num : vertices_list) {
//            INFO("here");
            size_t reverse_cost = size_t(-1);
            size_t intersected_cost = size_t(-1);
            size_t default_cost = tree.GetParentEdgeLength(clone_num);
//            INFO("here 2");
            if (best_reverse_edges.find(clone_num) != best_reverse_edges.end()) {
                reverse_cost = best_reverse_edges[clone_num]->Length();
            }
            if (best_intersected_edges.find(clone_num) != best_intersected_edges.end()) {
//                INFO("here 3");
                auto edge = best_intersected_edges[clone_num];
                if (SecondCloneIsFirstsAncestor(tree,
                                                edge->SrcNum(),
                                                clone_num)) {
                    continue;
                }
                if (edge->IsDoubleMutated()) {
                    intersected_cost = edge->Length();
                    if (intersected_cost < default_cost && intersected_cost < reverse_cost) {
                        tree.ReplaceEdge(clone_num, edge);
                        continue;
                    }
                } else {
//                    INFO("reconstruction case");
                    const auto &left = *edge->SrcClone();
                    const auto &right = *edge->DstClone();
                    size_t left_num = edge->SrcNum();
                    size_t right_num = clone_num;
                    size_t parent_num = clone_set.size();


                    auto cdr3_range = clone_by_read_constructor_.GetGeneCDR3Ranges(left.VGene(), left.JGene());
                    ++current_fake_clone_index_;
                    ++reconstructed_;
                    auto tpl = ParentReadReconstructor::ReconstructParentRead(left,
                                                                              right,
                                                                              current_fake_clone_index_,
                                                                              cdr3_range.first,
                                                                              cdr3_range.second);
                    auto parent_clone = clone_by_read_constructor_.GetCloneByReadAndAlignment(tpl,
                                                                                              left.VGene(),
                                                                                              left.JGene());
                    auto left_parent_edge = edge_constructor->ConstructEdge(left,
                                                                            parent_clone,
                                                                            left_num,
                                                                            parent_num);
                    auto parent_right_edge = edge_constructor->ConstructEdge(parent_clone,
                                                                             right,
                                                                             parent_num,
                                                                             right_num);
                    VERIFY(left_parent_edge->IsReverseDirected() && parent_right_edge->IsDirected());
                    VERIFY(left.CDR3Range().length() == right.CDR3Range().length() &&
                           left.CDR3Range().length() == parent_clone.CDR3Range().length());
                    intersected_cost = left_parent_edge->Length() + parent_right_edge->Length();
//                    INFO("reconstructed");
                    if (intersected_cost < default_cost && intersected_cost < reverse_cost) {
//                        bool transitive = true;
//                        clone_set.AddClone(parent_clone);
////                        INFO("inserting: intersected case");
//                        vertices_nums.insert(parent_num);
//                        tree.ReplaceEdge(parent_num, left_parent_edge);
//                        tree.ReplaceEdge(right_num, parent_right_edge);
//
//                        auto it = getRelatedClonesIterator(hamming_graph_info_, *edge->SrcClone());
//                        while (it.HasNext()) {
//                            size_t dst_num = it.Next();
//                            auto edge_n = edge_constructor->ConstructEdge(parent_clone,
//                                                                        clone_set[dst_num],
//                                                                        parent_num,
//                                                                        dst_num);
//                            if (edge_n->IsDirected() &&
//                                    (!tree.HasParentEdge(dst_num) ||
//                                     tree.GetParentEdgeLength(dst_num) > edge_n->Length()) &&
//                                    !SecondCloneIsFirstsAncestor(tree, parent_num, dst_num)) {
//                                tree.ReplaceEdge(dst_num, edge_n);
//                                INFO("replacing");
//                                transitive = false;
//                            }
//                        }
//                        continue;
                        tree.ReplaceEdge(clone_num, edge);
                    }
                    --current_fake_clone_index_;
                    --reconstructed_;
                }
            }
            if (reverse_cost < default_cost && reverse_cost < intersected_cost &&
                !SecondCloneIsFirstsAncestor(tree, best_reverse_edges[clone_num]->SrcNum(), clone_num)) {
//                INFO("reverse case");
                tree.ReplaceEdge(clone_num, best_reverse_edges[clone_num]);
                continue;
            }
        }
    }

    bool Base_CDR3_HG_CC_Processor::SecondCloneIsFirstsAncestor(EvolutionaryTree &tree,
                                                                size_t first_clone,
                                                                size_t second_clone) {
        size_t current_clone = first_clone;
        while (tree.HasParentEdge(current_clone)) {
            current_clone = tree.GetParentEdge(current_clone)->SrcNum();
            if (current_clone == second_clone) {
//                INFO("true");
                return true;
            }
//            INFO("mid");
        }
//        INFO("false");
        return false;
    }


    /*
    *       o
    *     /  \
    *    o    o
    *          \
    *           o <=  o
    *
    *
    *       ---->
    *
    *
    *           x
    *         /  \
    *       o     \
    *     /  \     \
    *    o    o     \
    *          \     \
    *           o     o
    */
    bool Base_CDR3_HG_CC_Processor::ReconstructAncestralLineageSimple(
            EvolutionaryEdgePtr edge,
            EvolutionaryTree &tree,
            boost::unordered_set<size_t> &vertices_nums,
            const std::shared_ptr<EvolutionaryEdgeConstructor> &edge_constructor,
            std::vector<size_t> &roots,
            boost::unordered_map<size_t, size_t> &iterator_index_map) {
        VERIFY_MSG(edge->IsIntersected(), "ancesrtal lineage reconstructor got a non-intersected edge");
        auto &clone_set = *clone_set_ptr_;
        size_t left_num = edge->DstNum();
        size_t right_num = edge->SrcNum();
        while (tree.HasParentEdge(left_num)) {
            left_num = tree.GetParentEdge(left_num)->SrcNum();
        }
        const auto &left = clone_set[left_num];
        const auto &right = clone_set[right_num];
        auto edge_n = edge_constructor->ConstructEdge(left, right, left_num, right_num);
        if (edge_n->IsDirected()) {
            tree.ReplaceEdge(right_num, edge_n);
            return false;
        }
        if (!CheckClonesConsistencyForReconstruction(left, right)) {
            return false;
        }
        size_t parent_num = clone_set.size();
        ++current_fake_clone_index_;
        ++reconstructed_;
        auto cdr3_range = clone_by_read_constructor_.GetGeneCDR3Ranges(left.VGene(), left.JGene());
        auto tpl = ParentReadReconstructor::ReconstructParentRead(left,
                                                                  right,
                                                                  current_fake_clone_index_,
                                                                  cdr3_range.first,
                                                                  cdr3_range.second);
        auto parent_clone = clone_by_read_constructor_.GetCloneByReadAndAlignment(tpl,
                                                                                  left.VGene(),
                                                                                  left.JGene());
        auto new_left_parent_edge = edge_constructor->ConstructEdge(parent_clone,
                                                                    left,
                                                                    parent_num,
                                                                    left_num);
        auto new_right_parent_edge = edge_constructor->ConstructEdge(parent_clone,
                                                                     right,
                                                                     parent_num,
                                                                     right_num);
        iterator_index_map[clone_set.size()] = left_num;

        VERIFY_MSG(new_left_parent_edge->IsDirected() && new_right_parent_edge->IsDirected(),
                   "error: edge from reconstructed parent is not directed");
        size_t left_cdr3_length = left.CDR3Range().length();
        size_t right_cdr3_length = right.CDR3Range().length();
        size_t parent_cdr3_length = parent_clone.CDR3Range().length();
        VERIFY(left_cdr3_length == right_cdr3_length && left_cdr3_length == parent_cdr3_length);
        clone_set.AddClone(parent_clone);
        vertices_nums.insert(parent_num);
        tree.ReplaceEdge(left_num, new_left_parent_edge);
        tree.ReplaceEdge(right_num, new_right_parent_edge);
        roots.push_back(parent_num);
        return true;
    }


    void Base_CDR3_HG_CC_Processor::HandleRootNeighbour(
            size_t root_num,
            size_t dst_num,
            boost::unordered_set<size_t> &vertices_nums,
            EvolutionaryTree &tree,
            boost::unordered_map<size_t, EvolutionaryEdgePtr> &roots_nearest_neighbours,
            const std::shared_ptr<EvolutionaryEdgeConstructor> &edge_constructor) {

        const auto &clone_set = *clone_set_ptr_;
        if (dst_num == root_num || vertices_nums.find(dst_num) == vertices_nums.end()) {
            return;
        }
        auto edge = edge_constructor->ConstructEdge(
                clone_set[root_num],
                clone_set[dst_num],
                root_num,
                dst_num);
        auto edge_r = edge_constructor->ConstructEdge(
                clone_set[dst_num],
                clone_set[root_num],
                dst_num,
                root_num);
        if (edge->IsDirected()) {
            if (!tree.HasParentEdge(dst_num) or edge->Length() < tree.GetParentEdge(dst_num)->Length()) {
                tree.ReplaceEdge(dst_num, edge);
            }
        }
        if (edge_r->IsDirected()) {
            if (!tree.HasParentEdge(root_num) or edge_r->Length() < tree.GetParentEdge(root_num)->Length()) {
                tree.ReplaceEdge(root_num, edge_r);
            }
        }
        if (edge->IsIntersected()) {
            if (roots_nearest_neighbours.find(root_num) == roots_nearest_neighbours.end() ||
                edge->Length() < roots_nearest_neighbours[root_num]->Length()) {
                roots_nearest_neighbours[root_num] = edge;
            }
        }
    }

    bool Base_CDR3_HG_CC_Processor::CheckClonesConsistencyForReconstruction(
            const annotation_utils::AnnotatedClone &left,
            const annotation_utils::AnnotatedClone &right) {
        if (left.CDR3Range().length() != right.CDR3Range().length()) {
            return false;
        }
        if (!annotation_utils::SHMComparator::SHMsInsertionBlocksAreEqual(left.VSHMs(),
                                                                          right.VSHMs()) ||
            !annotation_utils::SHMComparator::SHMsInsertionBlocksAreEqual(left.JSHMs(),
                                                                          right.JSHMs())) {
            return false;
        }
        return (left.CDR3Range().start_pos <= left.VAlignment().EndQueryPosition() + 1 &&
                left.CDR3Range().end_pos > left.VAlignment().EndQueryPosition() &&
                left.CDR3Range().start_pos <= left.JAlignment().StartQueryPosition() &&
                left.CDR3Range().end_pos >= left.JAlignment().StartQueryPosition() - 1 &&
                left.VAlignment().EndQueryPosition() < left.JAlignment().StartQueryPosition() &&
                right.CDR3Range().start_pos <= right.VAlignment().EndQueryPosition() + 1 &&
                right.CDR3Range().end_pos > right.VAlignment().EndQueryPosition() &&
                right.CDR3Range().start_pos <= right.JAlignment().StartQueryPosition() &&
                right.CDR3Range().end_pos >= right.JAlignment().StartQueryPosition() - 1 &&
                right.VAlignment().EndQueryPosition() < right.JAlignment().StartQueryPosition());
    }
}