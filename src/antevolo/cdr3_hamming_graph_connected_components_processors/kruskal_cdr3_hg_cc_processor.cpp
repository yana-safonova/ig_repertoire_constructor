#include <vj_class_processors/edmonds_tarjan_DMST_calculator.hpp>
#include "kruskal_cdr3_hg_cc_processor.hpp"
#include "parent_read_reconstructor.hpp"


namespace antevolo {
    void Kruskal_CDR3_HG_CC_Processor::SetUndirectedComponentsParentEdges(
            SparseGraphPtr hg_component,
            size_t component_id,
            boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges) {
        auto edge_constructor = GetEdgeConstructor();
        // adding directed edges between identical CDR3s
        for (size_t i = 0; i < hg_component->N(); i++) {
            size_t old_index = graph_component_.GetOldVertexByNewVertex(component_id, i);
            //auto clones_sharing_cdr3 = unique_cdr3s_map_[unique_cdr3s_[old_index]];
            auto clones_sharing_cdr3 = unique_cdr3s_map_.find(unique_cdr3s_[old_index])->second;
            for (size_t it1 = 0; it1 < clones_sharing_cdr3.size(); it1++)
                for(size_t it2 = it1 + 1; it2 < clones_sharing_cdr3.size(); it2++) {
                    size_t src_num = clones_sharing_cdr3[it1];
                    size_t dst_num = clones_sharing_cdr3[it2];
                    auto edge = edge_constructor->ConstructEdge(
                            clone_set_[src_num],
                            clone_set_[dst_num],
                            src_num,
                            dst_num);
                    SetUndirectedComponentParentEdge(ds_on_undirected_edges.find_set(dst_num),
                                                          edge);
                    auto edge_r = edge_constructor->ConstructEdge(
                            clone_set_[dst_num],
                            clone_set_[src_num],
                            dst_num,
                            src_num);
                    SetUndirectedComponentParentEdge(ds_on_undirected_edges.find_set(src_num),
                                                          edge_r);
                }
        }
        // adding directed edges between similar CDR3s
        for(size_t i = 0; i < hg_component->N(); i++) {
            size_t old_index1 = graph_component_.GetOldVertexByNewVertex(component_id, i);
            //for (size_t j = hg_component->RowIndex()[i]; j < hg_component->RowIndex()[i + 1]; j++) {
            for (auto it = hg_component->VertexEdges(i).begin(); it != hg_component->VertexEdges(i).end(); it++) {
                //size_t old_index2 = graph_component_.GetOldVertexByNewVertex(component_id, hg_component->Col()[j]);
                size_t old_index2 = graph_component_.GetOldVertexByNewVertex(component_id, *it);
                auto indices_1 = unique_cdr3s_map_.find(unique_cdr3s_[old_index1])->second;
                auto indices_2 = unique_cdr3s_map_.find(unique_cdr3s_[old_index2])->second;
                for (auto it1 = indices_1.begin(); it1 != indices_1.end(); it1++)
                    for (auto it2 = indices_2.begin(); it2 != indices_2.end(); it2++) {
                        auto edge = edge_constructor->ConstructEdge(
                                clone_set_[*it1],
                                clone_set_[*it2],
                                *it1,
                                *it2);
                        SetUndirectedComponentParentEdge(ds_on_undirected_edges.find_set(*it2),
                                                         edge);
//                        auto edge_r = edge_constructor->ConstructEdge(
//                                clone_set_[*it2],
//                                clone_set_[*it1],
//                                *it2,
//                                *it1);
//                        SetUndirectedComponentParentEdge(ds_on_undirected_edges.find_set(*it1),
//                                                              edge_r);
                    }
            }
        }
    }

    void Kruskal_CDR3_HG_CC_Processor::SetDirections(const boost::unordered_set<size_t>& vertices_nums,
                                                      EvolutionaryTree& tree,
                                                      boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges) {
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
                PrepareSubtreeKruskal(edge_vector, root, clone_set_, edge_constructor);
            }
            else {
                PrepareSubtreeKruskal(edge_vector, vertex.first, clone_set_, edge_constructor);
            }
            for (auto p : edge_vector) {
                auto edge = edge_constructor->ConstructEdge(
                        clone_set_[p.first],
                        clone_set_[p.second],
                        p.first,
                        p.second);
                tree.AddUndirected(p.second, edge);
            }
        }
    }

    void Kruskal_CDR3_HG_CC_Processor::ReconstructMissingVertices(const boost::unordered_set<size_t> &vertices_nums,
                                                                  EvolutionaryTree &tree,
                                                                  SparseGraphPtr hg_component,
                                                                  size_t component_id) {
        INFO("Reconstruction of missing vertices started");
        auto edge_constructor = GetEdgeConstructor();
        boost::unordered_map<size_t, EvolutionaryEdgePtr> roots_nearest_neighbours;
        std::vector<size_t> roots;
        for (auto v : vertices_nums) {
            if (!tree.Contains(v)) {
                roots.push_back(v);
            }
        }


        for (size_t r_i = 0; r_i < roots.size(); ++r_i) {
            INFO("Reconstructing parents for root " << r_i );
            size_t root_num = roots[r_i];
            std::string root_cdr3;
            const auto& root_cdr3_dna5 = clone_set_[root_num].CDR3();
            size_t  cdr3_length = seqan::length(root_cdr3_dna5);
            for (size_t i = 0; i < cdr3_length; ++i) {
                root_cdr3.push_back(root_cdr3_dna5[i]);
            }
            size_t old_index = cdr3_to_index_map_.find(root_cdr3)->second;
            // same cdr3
            const auto& clones_sharing_cdr3 =  unique_cdr3s_map_.find(root_cdr3)->second;
            HandleRoot(root_num,
                       unique_cdr3s_map_.find(root_cdr3)->second,
                       tree,
                       roots_nearest_neighbours,
                       edge_constructor);
            // similar cdr3
            size_t new_index = graph_component_.GetNewVertexByOldVertex(old_index);
            for (auto it = hg_component->VertexEdges(new_index).begin();
                 it != hg_component->VertexEdges(new_index).end();
                 it++) {

                size_t old_index2 = graph_component_.GetOldVertexByNewVertex(component_id, *it);
                HandleRoot(root_num,
                           unique_cdr3s_map_.find(unique_cdr3s_[old_index2])->second,
                           tree,
                           roots_nearest_neighbours,
                           edge_constructor);
            }
//            INFO("Reconstructing missing vertices before edge");
            if (roots_nearest_neighbours.find(root_num) == roots_nearest_neighbours.end()) {
                continue;
            }
            auto edge = roots_nearest_neighbours[root_num];
//            INFO("Reconstructing Ancestral Lineage from" <<
//                                                         edge->DstClone()->Read().name << "  " <<
//                                                         edge->SrcClone()->Read().name);
            ReconstructAncestralLineage(edge, tree, edge_constructor, roots);
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
                                                 CloneSetWithFakes& clone_set,
                                                 std::shared_ptr<EvolutionaryEdgeConstructor> edge_constructor) {
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

    void Kruskal_CDR3_HG_CC_Processor::ReconstructAncestralLineage(
            EvolutionaryEdgePtr edge,
            EvolutionaryTree& tree,
            const std::shared_ptr<EvolutionaryEdgeConstructor>& edge_constructor,
            std::vector<size_t>& roots) {
        while (tree.Contains(edge->SrcNum())) {
            VERIFY_MSG(edge->IsIntersected(), "ancesrtal lineage reconstructor got a non-intersected edge");
            const auto& left = edge->SrcClone();
            const auto& right = edge->DstClone();
            INFO("creating cparent from" << left->Read().name << "  " << right->Read().name);
            auto parent_read = ParentReadReconstructor::ReconstructParentRead(left, right, clone_set_.size());
            auto parent_clone = clone_by_read_constructor_.GetCloneByRead(parent_read);
            auto old_left_parent_edge = tree.GetParentEdge(edge->SrcNum());
            auto new_left_parent_edge = edge_constructor->ConstructEdge(parent_clone,
                                                                        *left,
                                                                        clone_set_.size(),
                                                                        edge->SrcNum());
            if (new_left_parent_edge->Length() > old_left_parent_edge->Length()) {
                break;
            }
            auto new_right_parent_edge = edge_constructor->ConstructEdge(parent_clone,
                                                                         *right,
                                                                         clone_set_.size(),
                                                                         edge->SrcNum());
            VERIFY_MSG(new_left_parent_edge->IsDirected() && new_right_parent_edge->IsDirected(),
                       "error: edge from reconstructed parent is not directed");
            clone_set_.AddClone(parent_clone);
            edge = edge_constructor->ConstructEdge(*old_left_parent_edge->SrcClone(),
                                                   clone_set_[clone_set_.size()],
                                                   old_left_parent_edge->SrcNum(),
                                                   clone_set_.size());
            tree.ReplaceEdge(old_left_parent_edge->DstNum(), new_left_parent_edge);
            tree.ReplaceEdge(clone_set_.size(), new_right_parent_edge);
        }

        size_t left_num = edge->SrcNum();
        size_t right_num = edge->DstNum();
        while (left_num) {
            left_num = tree.GetParentEdge(left_num)->SrcNum();
        }
        edge = edge_constructor->ConstructEdge(clone_set_[left_num],
                                               clone_set_[right_num],
                                               left_num,
                                               right_num);
        const auto& left = edge->SrcClone();
        const auto& right = edge->DstClone();
        auto parent_read = ParentReadReconstructor::ReconstructParentRead(left, right, clone_set_.size());
        auto parent_clone = clone_by_read_constructor_.GetCloneByRead(parent_read);
        auto new_left_parent_edge = edge_constructor->ConstructEdge(parent_clone,
                                                                    *left,
                                                                    clone_set_.size(),
                                                                    left_num);
        auto new_right_parent_edge = edge_constructor->ConstructEdge(parent_clone,
                                                                     *right,
                                                                     clone_set_.size(),
                                                                     right_num);
        VERIFY_MSG(new_left_parent_edge->IsDirected() && new_right_parent_edge->IsDirected(),
                   "error: edge from reconstructed parent is not directed");
        clone_set_.AddClone(parent_clone);
        tree.ReplaceEdge(left_num, new_left_parent_edge);
        tree.ReplaceEdge(right_num, new_right_parent_edge);
        roots.push_back(clone_set_.size());
    }

    void Kruskal_CDR3_HG_CC_Processor::HandleRoot(
            size_t root_num,
            const std::vector<size_t>& clones_sharing_cdr3,
            EvolutionaryTree& tree,
            boost::unordered_map<size_t, EvolutionaryEdgePtr>& roots_nearest_neighbours,
            const std::shared_ptr<EvolutionaryEdgeConstructor>& edge_constructor){

        for (size_t i = 0; i < clones_sharing_cdr3.size(); ++i) {
//            INFO("handling root " << root_num << " it " << i);
            size_t neighbour_num = clones_sharing_cdr3[i];
            if (root_num == neighbour_num) {
                continue;
            }
            if (clone_set_.IsFake(root_num)) {
                auto edge_r = edge_constructor->ConstructEdge(
                        clone_set_[root_num],
                        clone_set_[neighbour_num],
                        root_num,
                        neighbour_num);
                if (edge_r->IsDirected()) {
                    tree.AddDirected(neighbour_num, edge_r);
                }
            }
            auto edge = edge_constructor->ConstructEdge(
                    clone_set_[neighbour_num],
                    clone_set_[root_num],
                    neighbour_num,
                    root_num);
            VERIFY_MSG(edge->IsIntersected(), "missing vertices reconstructor got a non-intersected edge");
            if (roots_nearest_neighbours.find(root_num) == roots_nearest_neighbours.end() ||
                roots_nearest_neighbours[root_num]->Length() > edge->Length()) {
                roots_nearest_neighbours[root_num] = edge;
//                INFO("found a better edge for root " << root_num);
            }
        }
//        INFO("end handling root " << root_num);
    }

}