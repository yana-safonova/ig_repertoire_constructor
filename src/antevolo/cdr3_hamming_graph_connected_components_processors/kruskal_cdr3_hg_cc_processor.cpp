#include <vj_class_processors/edmonds_tarjan_DMST_calculator.hpp>
#include "kruskal_cdr3_hg_cc_processor.hpp"
#include "parent_read_reconstructor.hpp"


namespace antevolo {
    const size_t EVO_EDGE_MAX_LENGTH = 400;

    void Kruskal_CDR3_HG_CC_Processor::SetUndirectedComponentsParentEdges(
            boost::disjoint_sets<AP_map, AP_map>& ds_on_undirected_edges,
            const boost::unordered_set<size_t>& vertices_nums) {
        const auto& clone_set = *clone_set_ptr_;
        auto edge_constructor = GetEdgeConstructor();
        for (auto src_num : vertices_nums) {
            size_t dst_num;
            auto it = getRelatedClonesIterator(hamming_graph_info_, clone_set[src_num]);
            while (it.HasNext()) {
                dst_num = it.Next();
                VERIFY(vertices_nums.find(dst_num) != vertices_nums.end());
                if (dst_num == src_num) {
                    continue;
                }
                auto edge = edge_constructor->ConstructEdge(
                        clone_set[src_num],
                        clone_set[dst_num],
                        src_num,
                        dst_num);
                SetUndirectedComponentParentEdge(ds_on_undirected_edges.find_set(dst_num),
                                                 edge);
            }
        }
    }

    void Kruskal_CDR3_HG_CC_Processor::SetDirections(boost::disjoint_sets<AP_map, AP_map>& ds_on_undirected_edges,
                                                     const boost::unordered_set<size_t>& vertices_nums,
                                                     EvolutionaryTree &tree) {
        const auto& clone_set = *clone_set_ptr_;
        auto edge_constructor = GetEdgeConstructor();
        boost::unordered_set<size_t> undirected_graph_vertices;
        for (auto p : undirected_graph_) {
            undirected_graph_vertices.insert(p.first);
        }

        for (auto clone_num : vertices_nums) {
            if (undirected_graph_vertices.find(clone_num) == undirected_graph_vertices.end()) {
                // if it is an undirected-isolated vertex
                if (GetUndirectedCompopentRoot(ds_on_undirected_edges.find_set(clone_num)) != size_t(-1)) {
                    const EvolutionaryEdgePtr& edge = GetUndirectedComponentParentEdge(clone_num);
                    tree.AddDirected(clone_num, edge/*, model_*/);
                };

                continue;
            }
            if (parent_edge_handled_[clone_num]) {
                continue;
            }
            auto vertex_it = undirected_graph_.find(clone_num);
            auto vertex = *vertex_it;
            std::vector<std::pair<size_t, size_t>> edge_vector;
            size_t root = GetUndirectedCompopentRoot(ds_on_undirected_edges.find_set(vertex.first));
            if (root != size_t(-1)) {
                const EvolutionaryEdgePtr& edge = GetUndirectedComponentParentEdge(ds_on_undirected_edges.find_set(vertex.first));
                tree.AddDirected(edge->DstNum(), edge/*, model_*/);
                PrepareSubtreeKruskal(edge_vector, root, edge_constructor);
            }
            else {
                PrepareSubtreeKruskal(edge_vector, vertex.first, edge_constructor);
            }
            for (auto p : edge_vector) {
                auto edge = edge_constructor->ConstructEdge(
                        clone_set[p.first],
                        clone_set[p.second],
                        p.first,
                        p.second);
                tree.AddUndirected(p.second, edge);
            }
        }
    }

    void Kruskal_CDR3_HG_CC_Processor::ReconstructMissingVertices(boost::unordered_set<size_t>& vertices_nums,
                                                                  EvolutionaryTree& tree) {
        const auto& clone_set = *clone_set_ptr_;
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
            auto it = getRelatedClonesIterator(hamming_graph_info_, clone_set[root_num]);
            while (it.HasNext()) {
                size_t dst_num = it.Next();
                VERIFY(vertices_nums.find(dst_num) != vertices_nums.end());
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
            for (size_t i = 0; i < roots.size(); ++i)  {
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
            const auto& root_clone = clone_set[root_num];
            const auto& neighbour_clone = *edge->DstClone();

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

    void Kruskal_CDR3_HG_CC_Processor::PrepareSubtree(std::vector<std::pair<size_t, size_t>>& edge_vector,
                                                          size_t root_num) {
        if (undirected_graph_.find(root_num) != undirected_graph_.end() &&
                !parent_edge_handled_[root_num]) {
            parent_edge_handled_[root_num] = true;
            for (size_t u : undirected_graph_[root_num]) {
                if (!parent_edge_handled_[u]) {
                    edge_vector.push_back(std::make_pair(root_num, u));
                    PrepareSubtree(edge_vector, u);
                }
            }
        }
    }

    void Kruskal_CDR3_HG_CC_Processor::PrepareSubtreeVertices(
            boost::unordered_set<size_t>& vertices_set,
            size_t root_num) {
        //we assume that component size is > 1
        if (vertices_set.find(root_num) != vertices_set.end()) {
            return;
        }
        vertices_set.insert(root_num);
        for (size_t u : undirected_graph_[root_num]) {
            PrepareSubtreeVertices(vertices_set, u);
        }
    }

    void Kruskal_CDR3_HG_CC_Processor::PrepareSubtreeKruskal(std::vector<std::pair<size_t, size_t>>& edge_vector,
                                                 size_t root_vertex,
                                                 std::shared_ptr<EvolutionaryEdgeConstructor> edge_constructor) {
        const auto& clone_set = *clone_set_ptr_;
        boost::unordered_set<size_t> vertices_set;
        PrepareSubtreeVertices(vertices_set, root_vertex);
        for (size_t v : vertices_set) {
            undirected_graph_[v].clear();
        }
        size_t n = vertices_set.size();
        boost::unordered_map<size_t, size_t> vertex_to_index;
        std::vector<size_t> index_to_vertex(n);
        size_t index = 0;
        for (size_t vertex : vertices_set) {
            vertex_to_index[vertex] = index;
            index_to_vertex[index] = vertex;
            ++index;
        }
        typedef EdmondsTarjanDMSTCalculator::WeightedEdge WeightedEdge;
        std::vector<WeightedEdge> edges;
        std::vector<boost::unordered_map<size_t, size_t>> kruskal_graph(n);
        for (size_t v = 0; v < n; ++v) {
            for (size_t u = 0; u < n; ++u) {
                if (v == u) {
                    continue;
                }
                size_t CDR3_dist = edge_constructor->ConstructEdge(clone_set[index_to_vertex[v]],
                                                                   clone_set[index_to_vertex[u]],
                                                                   index_to_vertex[v],
                                                                   index_to_vertex[u])->CDR3Distance();
                if (CDR3_dist <= config_.similar_cdr3s_params.num_mismatches) {
                    edges.push_back(WeightedEdge(v, u, static_cast<double>(CDR3_dist)));
                }
            }
        }

        std::map<size_t, size_t> rank;
        std::map<size_t, size_t> parent;
        typedef boost::associative_property_map<std::map<size_t, size_t>> AP_map;
        boost::disjoint_sets<AP_map, AP_map> ds(
                boost::make_assoc_property_map(rank),
                boost::make_assoc_property_map(parent));
        for (size_t i = 0; i < n; ++i) {
            ds.make_set(i);
        }
        std::sort(edges.begin(), edges.end(), [](WeightedEdge e1, WeightedEdge e2) {
            return e1.weight_ < e2.weight_;
        });

        std::vector<WeightedEdge> edges_to_add;
        size_t edge_num = 0;
        size_t added_edge_num = 0;
        while (edge_num < edges.size() && added_edge_num < n) {
            auto const &edge = edges[edge_num];
            if (ds.find_set(edge.src_) == ds.find_set(edge.dst_)) {
                ++edge_num;
                continue;
            }
            edges_to_add.push_back(edge);
            ds.union_set(edge.src_, edge.dst_);
            ++edge_num;
            ++added_edge_num;
        }

        for (auto &we : edges_to_add) {
            size_t src_clone_num = index_to_vertex[we.src_];
            size_t dst_clone_num = index_to_vertex[we.dst_];
            undirected_graph_[src_clone_num].insert(dst_clone_num);
            undirected_graph_[dst_clone_num].insert(src_clone_num);
        }
        PrepareSubtree(edge_vector, root_vertex);
    }

    void Kruskal_CDR3_HG_CC_Processor::SetUndirectedComponentParentEdge(size_t root_num,
                                                                        EvolutionaryEdgePtr edge) {
        if(edge->IsDirected()) {
            if (undirected_components_edges_.find(root_num) == undirected_components_edges_.end()) {
                undirected_components_edges_[root_num] = edge;
                return;
            }
            const EvolutionaryEdgePtr& parent_edge =  undirected_components_edges_[root_num];
            if (parent_edge->Length() > edge->Length()) { // todo: compare only by added shms and then by cdr3?
                //if clone_set_[*it2] is root or if the new edge is shorter
                undirected_components_edges_[root_num] = edge;
                return;
            }
            /*
            if (parent_edge.num_added_shms == edge.num_added_shms && parent_edge.cdr3_distance > edge.cdr3_distance) {
                //if clone_set_[*it2] is root or if the new edge is shorter
                undirected_components_edges_[root_num] = edge;
                return;
            }
            */
        }
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
    bool Kruskal_CDR3_HG_CC_Processor::ReconstructAncestralLineageSimple(
            EvolutionaryEdgePtr edge,
            EvolutionaryTree &tree,
            boost::unordered_set<size_t> &vertices_nums,
            const std::shared_ptr<EvolutionaryEdgeConstructor> &edge_constructor,
            std::vector<size_t> &roots,
            boost::unordered_map<size_t, size_t> &iterator_index_map) {
        VERIFY_MSG(edge->IsIntersected(), "ancesrtal lineage reconstructor got a non-intersected edge");
        auto& clone_set = *clone_set_ptr_;
        size_t left_num = edge->DstNum();
        size_t right_num = edge->SrcNum();
        while (tree.HasParentEdge(left_num)) {
            left_num = tree.GetParentEdge(left_num)->SrcNum();
        }
        const auto& left = clone_set[left_num];
        const auto& right = clone_set[right_num];
        auto edge_n = edge_constructor->ConstructEdge(left, right, left_num, right_num);
        VERIFY(left.CDR3Range().length() == right.CDR3Range().length());
        if (edge_n->IsDirected()) {
            tree.ReplaceEdge(right_num, edge_n);
            return false;
        }
        if (!annotation_utils::SHMComparator::SHMsInsertionBlocksAreEqual(left.VSHMs(),
                                                                          right.VSHMs()) ||
            !annotation_utils::SHMComparator::SHMsInsertionBlocksAreEqual(left.JSHMs(),
                                                                          right.JSHMs())) {
            ++rejected_;
            return false;
        }
        if (left.CDR3Range().end_pos < left.VAlignment().EndQueryPosition() ||
            left.CDR3Range().start_pos > left.JAlignment().StartQueryPosition() ||
            right.CDR3Range().end_pos < right.VAlignment().EndQueryPosition() ||
            right.CDR3Range().start_pos > right.JAlignment().StartQueryPosition()) {
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
        auto& parent_read = parent_clone.Read();

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
     *       o     x
     *     /     /  \
     *    o    o     x
     *             /  \
     *           o     o
     */

//    void Kruskal_CDR3_HG_CC_Processor::ReconstructAncestralLineage(
//            EvolutionaryEdgePtr edge,
//            EvolutionaryTree& tree,
//            const std::shared_ptr<EvolutionaryEdgeConstructor>& edge_constructor,
//            std::vector<size_t>& roots) {
//        \}

    void Kruskal_CDR3_HG_CC_Processor::HandleRootNeighbour(
            size_t root_num,
            size_t dst_num,
            boost::unordered_set<size_t>& vertices_nums,
            EvolutionaryTree& tree,
            boost::unordered_map<size_t, EvolutionaryEdgePtr>& roots_nearest_neighbours,
            const std::shared_ptr<EvolutionaryEdgeConstructor>& edge_constructor) {

        const auto& clone_set = *clone_set_ptr_;
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


}