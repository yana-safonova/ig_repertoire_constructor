#include <vj_class_processors/edmonds_tarjan_DMST_calculator.hpp>
#include "kruskal_cdr3_hg_cc_processor.hpp"
#include "parent_read_reconstructor.hpp"


namespace antevolo {

    EvolutionaryTree Kruskal_CDR3_HG_CC_Processor::ConstructForest() {
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
        Refine(vertices_nums, tree);
        tree.AddAllEdges();
        return tree;
    }

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

}