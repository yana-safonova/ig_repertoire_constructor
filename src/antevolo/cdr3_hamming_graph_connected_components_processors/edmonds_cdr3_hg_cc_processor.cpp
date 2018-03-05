#include "edmonds_cdr3_hg_cc_processor.hpp"

namespace antevolo {
    EvolutionaryTree Edmonds_CDR3_HG_CC_Processor::Process() {
        EvolutionaryTree tree(clone_set_ptr_);
        vertices_nums_ = boost::unordered_set<size_t>(hamming_graph_info_.GetAllClones());

        for (auto it = vertices_nums_.begin(); it != vertices_nums_.end(); ) {
            if (clone_set_ptr_->operator[](*it).AnnotationIsNotValid()) {
                vertices_nums_.erase(it++);
            }
            else {
                ++it;
            }
        }

//        size_t cdr3_length = clone_set_ptr_->operator[](*vertices_nums_.cbegin()).CDR3Range().length();
//        for (size_t clone_num : vertices_nums_) {
//            VERIFY(clone_set_ptr_->operator[](clone_num).CDR3Range().length() == cdr3_length);
//        }

        for (size_t clone_num : vertices_nums_) {
//            VERIFY(clone_set_ptr_->operator[](clone_num).CDR3Range().length() == cdr3_length);
            tree.AddVertex(clone_num);
        }

        SetShortestDirectedParentEdges();
        auto input_edges = PrepareEdgeVector();
        auto branching_edges = EdmondsProcessor().process_edge_list(input_edges);
        SetEdges(tree, branching_edges);
        ReconstructMissingVertices(vertices_nums_, tree);
        Refine(vertices_nums_, tree);
        tree.AddAllEdges();
        return tree;
    }

    void Edmonds_CDR3_HG_CC_Processor::SetShortestDirectedParentEdges() {
        const auto &clone_set = *clone_set_ptr_;
        auto edge_constructor = GetEdgeConstructor();

        for (auto src_num : vertices_nums_) {
//            size_t dst_num;
//            auto it = getRelatedClonesIterator(hamming_graph_info_, clone_set[src_num]);
//            while (it.HasNext()) {
//                dst_num = it.Next();
//                VERIFY(vertices_nums_.find(dst_num) != vertices_nums_.end());

            for (size_t dst_num : vertices_nums_) {
                if (dst_num == src_num) {
                    continue;
                }
                auto edge = edge_constructor->ConstructEdge(clone_set[src_num],
                                                            clone_set[dst_num],
                                                            src_num,
                                                            dst_num);
                if (!edge->IsDirected()) {
                    continue;
                }
                if (shorthest_directed_edge_.find(dst_num) == shorthest_directed_edge_.end() ||
                    //                    GetLength(edge) < GetLength(shorthest_directed_edge_[dst_num])) {
                    edge->Length() < shorthest_directed_edge_[dst_num]->Length()) {
                    shorthest_directed_edge_[dst_num] = edge;
                }
            }
        }
    }

    std::vector<WeightedEdge<int>> Edmonds_CDR3_HG_CC_Processor::PrepareEdgeVector() {
        const auto &clone_set = *clone_set_ptr_;
        auto edge_constructor = GetEdgeConstructor();

        size_t germline_vertex = size_t(-1);

        std::vector<WeightedEdge<int>> res;

        for (auto src_num : vertices_nums_) {
            const auto &src_clone = clone_set[src_num];
            res.push_back(WeightedEdge<int>(germline_vertex, src_num,
                                            static_cast<int>(src_clone.VSHMs().size() + src_clone.JSHMs().size())));
//            size_t dst_num;
//            auto it = getRelatedClonesIterator(hamming_graph_info_, clone_set[src_num]);
//            while (it.HasNext()) {
//                dst_num = it.Next();
//                VERIFY(vertices_nums_.find(dst_num) != vertices_nums_.end());
            for (size_t dst_num : vertices_nums_) {
                if (dst_num == src_num) {
                    continue;
                }
                auto edge = edge_constructor->ConstructEdge(clone_set[src_num],
                                                            clone_set[dst_num],
                                                            src_num,
                                                            dst_num);
                if (edge->IsUndirected()) {
//                    res.push_back(EdmondsProcessor::WeightedEdge(src_num, dst_num, GetLength(edge)));
                    res.push_back(WeightedEdge<int>(src_num, dst_num, static_cast<int>(edge->Length())));
                }
            }
        }
        for (auto p : shorthest_directed_edge_) {
            auto edge = p.second;
//            res.push_back(EdmondsProcessor::WeightedEdge(edge->SrcNum(), edge->DstNum(), GetLength(edge)));
            res.push_back(WeightedEdge<int>(edge->SrcNum(), edge->DstNum(), static_cast<int>(edge->Length())));
        }
        return res;
    }

    void Edmonds_CDR3_HG_CC_Processor::SetEdges(EvolutionaryTree &tree,
                                                const std::vector<WeightedEdge<int>> &edge_vector) {
        const auto &clone_set = *clone_set_ptr_;
        auto edge_constructor = GetEdgeConstructor();

        for (auto we : edge_vector) {
//            std::cout << we.src_ << " " << we.dst_ << "\n";
            if (we.src_ == size_t(-1)) {
                continue;
            }
            auto edge = edge_constructor->ConstructEdge(clone_set[we.src_],
                                                        clone_set[we.dst_],
                                                        we.src_,
                                                        we.dst_);
            VERIFY(edge->IsDirected() || edge->IsUndirected());
            tree.ReplaceEdge(we.dst_, edge);
        }
    }
}