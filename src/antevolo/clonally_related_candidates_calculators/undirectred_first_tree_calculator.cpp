#include <path_helper.hpp>
#include "undirectred_first_tree_calculator.hpp"
#include "../../graph_utils/graph_io.hpp"
#include "../../graph_utils/graph_splitter.hpp"
#include <convert.hpp>
#include <annotation_utils/shm_comparator.hpp>
#include <evolutionary_graph_utils/evolutionary_graph_constructor.hpp>



namespace antevolo {
    void UndirectedFirstTreeCalculator::Clear() {
        unique_cdr3s_.clear();
        unique_cdr3s_map_.clear();
    }

    void UndirectedFirstTreeCalculator::CreateUniqueCDR3Map(
            core::DecompositionClass decomposition_class) {
        for(auto it = decomposition_class.begin(); it != decomposition_class.end(); it++) {
            if(clone_set_[*it].RegionIsEmpty(annotation_utils::StructuralRegion::CDR3))
                continue;
            auto cdr3 = core::dna5String_to_string(clone_set_[*it].CDR3());
            if(unique_cdr3s_map_.find(cdr3) == unique_cdr3s_map_.end())
                unique_cdr3s_map_[cdr3] = std::vector<size_t>();
            unique_cdr3s_map_[cdr3].push_back(*it);
        }
        for(auto it = unique_cdr3s_map_.begin(); it != unique_cdr3s_map_.end(); it++)
            unique_cdr3s_.push_back(it->first);
    }

    std::string UndirectedFirstTreeCalculator::GetFastaFname(core::DecompositionClass decomposition_class) {
        std::stringstream ss;
        size_t key = *decomposition_class.begin();
        ss << "CDR3_sequences_key_" << key << ".fasta";
        return path::append_path(output_params_.cdr_graph_dir, ss.str());
    }

    std::string UndirectedFirstTreeCalculator::GetGraphFname(core::DecompositionClass decomposition_class) {
        std::stringstream ss;
        size_t key = *decomposition_class.begin();
        ss << "CDR3_sequences_key_" << key << ".graph";
        return path::append_path(output_params_.cdr_graph_dir, ss.str());
    }

    std::string UndirectedFirstTreeCalculator::WriteUniqueCDR3InFasta(core::DecompositionClass decomposition_class) {
        std::string output_fname = GetFastaFname(decomposition_class);
        std::ofstream out(output_fname);
        for(auto it = unique_cdr3s_.begin(); it != unique_cdr3s_.end(); it++)
            out << ">" << *it << std::endl << *it << std::endl;
        return output_fname;
    }

    // return connected components of Hamming graph on CDR3s
    std::vector<SparseGraphPtr> UndirectedFirstTreeCalculator::ComputeCDR3HammingGraphs(std::string cdr_fasta,
                                                                                         std::string graph_fname) {
        std::string run_graph_constructor = "./build/release/bin/ig_swgraph_construct";
        std::stringstream ss;
        ss << run_graph_constructor << " -i " << cdr_fasta <<
                " -o " << graph_fname << " --tau " << num_mismatches_ << " -S " << " 0 " <<
                " -T " << " 0 " << " -k 10 > " << output_params_.trash_output;
        int err_code = system(ss.str().c_str());
        VERIFY_MSG(err_code == 0, "Graph constructor finished abnormally, error code: " << err_code);
        auto sparse_cdr_graph_ = GraphReader(graph_fname).CreateGraph();
        TRACE("Hamming graph contains " << sparse_cdr_graph_->N() << " edges and " << sparse_cdr_graph_->NZ() << " edges");
        //if(sparse_cdr_graph_->NZ() != 0)
        //    std::cout << *sparse_cdr_graph_ << std::endl;
        auto connected_components = ConnectedComponentGraphSplitter(sparse_cdr_graph_).Split();
        graph_component_ = sparse_cdr_graph_->GetGraphComponentMap();
        return connected_components;
    }

    void UndirectedFirstTreeCalculator::AddComponent(SparseGraphPtr hg_component,
                                                     size_t component_id, EvolutionaryTree& tree) {
        boost::unordered_set<size_t> vertices_nums;
        for (size_t i = 0; i < hg_component->N(); i++) {
            size_t old_index = graph_component_.GetOldVertexByNewVertex(component_id, i);
            auto clones_sharing_cdr3 = unique_cdr3s_map_[unique_cdr3s_[old_index]];
            for (size_t clone_num : clones_sharing_cdr3) {
                vertices_nums.insert(clone_num);
            }
        }

        //typedef std::map<size_t, size_t> AP_map;
        std::map<size_t, size_t> rank;
        std::map<size_t, size_t> parent;
        boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges(rank, parent);
        for (size_t i : vertices_nums) {
            ds_on_undirected_edges.make_set(i);
        }
        AddUndirectedForestToTheTree(hg_component, component_id, tree, ds_on_undirected_edges);
        //AddComponentToTheTree(hg_component, component_id, tree);
        SetUndirectedComponentsParentEdges(hg_component, component_id, tree, ds_on_undirected_edges);
        SetDirections(vertices_nums, tree, ds_on_undirected_edges);
    }

    void UndirectedFirstTreeCalculator::AddUndirectedForestToTheTree(SparseGraphPtr hg_component,
                                                                     size_t component_id, EvolutionaryTree& tree,
                                                                     boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges) {

        //adding undirected edges first
        for (size_t i = 0; i < hg_component->N(); i++) {
            size_t old_index = graph_component_.GetOldVertexByNewVertex(component_id, i);
            auto clones_sharing_cdr3 = unique_cdr3s_map_[unique_cdr3s_[old_index]];
            for (size_t it1 = 0; it1 < clones_sharing_cdr3.size(); it1++)
                for (size_t it2 = it1 + 1; it2 < clones_sharing_cdr3.size(); it2++) {
                    size_t src_num = clones_sharing_cdr3[it1];
                    size_t dst_num = clones_sharing_cdr3[it2];
                    // if clones are not in the same connected component and
                    // the edge is undirected
                    if (ds_on_undirected_edges.find_set(src_num) !=
                        ds_on_undirected_edges.find_set(dst_num) &&
                        annotation_utils::SHMComparator::SHMsAreEqual(
                                clone_set_[src_num].VSHMs(),
                                clone_set_[dst_num].VSHMs())) {
                        tree.AddUndirectedPair(src_num, dst_num);
                        ds_on_undirected_edges.union_set(src_num, dst_num);
                    }
                }
        }
        for (size_t i = 0; i < hg_component->N(); i++)
            for (size_t j = hg_component->RowIndex()[i]; j < hg_component->RowIndex()[i + 1]; j++) {
                size_t old_index1 = graph_component_.GetOldVertexByNewVertex(component_id, i);
                size_t old_index2 = graph_component_.GetOldVertexByNewVertex(component_id, hg_component->Col()[j]);
                auto indices_1 = unique_cdr3s_map_[unique_cdr3s_[old_index1]];
                auto indices_2 = unique_cdr3s_map_[unique_cdr3s_[old_index2]];
                for (auto it1 = indices_1.begin(); it1 != indices_1.end(); it1++)
                    for (auto it2 = indices_2.begin(); it2 != indices_2.end(); it2++) {
                        if (ds_on_undirected_edges.find_set(*it1) !=
                            ds_on_undirected_edges.find_set(*it2) &&
                            annotation_utils::SHMComparator::SHMsAreEqual(
                                    clone_set_[*it1].VSHMs(),
                                    clone_set_[*it2].VSHMs())) {
                            tree.AddUndirectedPair(*it1, *it2);
                            ds_on_undirected_edges.union_set(*it1, *it2);
                        }
                    }
            }
    }
    void UndirectedFirstTreeCalculator::AddComponentToTheTree(SparseGraphPtr hg_component,
                                                                     size_t component_id, EvolutionaryTree& tree) {
        auto edge_constructor = GetEdgeConstructor();
        // adding directed edges between identical CDR3s
        for(size_t i = 0; i < hg_component->N(); i++) {
            size_t old_index = graph_component_.GetOldVertexByNewVertex(component_id, i);
            auto clones_sharing_cdr3 = unique_cdr3s_map_[unique_cdr3s_[old_index]];
            for(size_t it1 = 0; it1 < clones_sharing_cdr3.size(); it1++)
                for(size_t it2 = it1 + 1; it2 < clones_sharing_cdr3.size(); it2++) {
                    size_t src_num = clones_sharing_cdr3[it1];
                    size_t dst_num = clones_sharing_cdr3[it2];
                    auto edge = edge_constructor->ConstructEdge(
                            clone_set_[src_num],
                            clone_set_[dst_num],
                            src_num,
                            dst_num);
                    tree.AddDirected(dst_num, edge, model_);
                    std::vector<std::pair<size_t, size_t>> edge_vector;
                    tree.PrepareSubtree(edge_vector, dst_num);
                    for (auto p : edge_vector) {
                        auto edge2 = edge_constructor->ConstructEdge(
                                clone_set_[p.first],
                                clone_set_[p.second],
                                p.first,
                                p.second);
                        tree.AddUndirected(p.second, edge2);
                    }
                }
        }
        // adding directed edges between similar CDR3s
        for(size_t i = 0; i < hg_component->N(); i++)
            for(size_t j = hg_component->RowIndex()[i]; j < hg_component->RowIndex()[i + 1]; j++) {
                size_t old_index1 = graph_component_.GetOldVertexByNewVertex(component_id, i);
                size_t old_index2 = graph_component_.GetOldVertexByNewVertex(component_id, hg_component->Col()[j]);
                auto indices_1 = unique_cdr3s_map_[unique_cdr3s_[old_index1]];
                auto indices_2 = unique_cdr3s_map_[unique_cdr3s_[old_index2]];
                for(auto it1 = indices_1.begin(); it1!= indices_1.end(); it1++)
                    for(auto it2 = indices_2.begin(); it2!= indices_2.end(); it2++) {
                        auto edge = edge_constructor->ConstructEdge(
                                clone_set_[*it1],
                                clone_set_[*it2],
                                *it1,
                                *it2);
                        tree.AddDirected(*it2, edge, model_);
                        std::vector<std::pair<size_t, size_t>> edge_vector;
                        tree.PrepareSubtree(edge_vector, *it2);
                        for (auto p : edge_vector) {
                            auto edge2 = edge_constructor->ConstructEdge(
                                    clone_set_[p.first],
                                    clone_set_[p.second],
                                    p.first,
                                    p.second);
                            tree.AddUndirected(p.second, edge2);
                        }
                    }
            }
        auto undirected_graph = tree.GetUndirectedGraph();
        for (auto edges : undirected_graph) {
            std::vector<std::pair<size_t, size_t>> edge_vector;
            tree.PrepareSubtree(edge_vector, edges.first);
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

    void UndirectedFirstTreeCalculator::SetUndirectedComponentsParentEdges(SparseGraphPtr hg_component,
                                                                           size_t component_id, EvolutionaryTree& tree,
                                                                           boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges) {
        auto edge_constructor = GetEdgeConstructor();
        // adding directed edges between identical CDR3s
        for(size_t i = 0; i < hg_component->N(); i++) {
            size_t old_index = graph_component_.GetOldVertexByNewVertex(component_id, i);
            auto clones_sharing_cdr3 = unique_cdr3s_map_[unique_cdr3s_[old_index]];
            for(size_t it1 = 0; it1 < clones_sharing_cdr3.size(); it1++)
                for(size_t it2 = it1 + 1; it2 < clones_sharing_cdr3.size(); it2++) {
                    size_t src_num = clones_sharing_cdr3[it1];
                    size_t dst_num = clones_sharing_cdr3[it2];
                    auto edge = edge_constructor->ConstructEdge(
                            clone_set_[src_num],
                            clone_set_[dst_num],
                            src_num,
                            dst_num);
                    tree.SetUndirectedComponentParentEdge(ds_on_undirected_edges.find_set(dst_num),
                                                          edge, model_);
                }
        }
        // adding directed edges between similar CDR3s
        for(size_t i = 0; i < hg_component->N(); i++)
            for(size_t j = hg_component->RowIndex()[i]; j < hg_component->RowIndex()[i + 1]; j++) {
                size_t old_index1 = graph_component_.GetOldVertexByNewVertex(component_id, i);
                size_t old_index2 = graph_component_.GetOldVertexByNewVertex(component_id, hg_component->Col()[j]);
                auto indices_1 = unique_cdr3s_map_[unique_cdr3s_[old_index1]];
                auto indices_2 = unique_cdr3s_map_[unique_cdr3s_[old_index2]];
                for(auto it1 = indices_1.begin(); it1!= indices_1.end(); it1++)
                    for(auto it2 = indices_2.begin(); it2!= indices_2.end(); it2++) {
                        auto edge = edge_constructor->ConstructEdge(
                                clone_set_[*it1],
                                clone_set_[*it2],
                                *it1,
                                *it2);
                        tree.SetUndirectedComponentParentEdge(ds_on_undirected_edges.find_set(*it2),
                                                              edge, model_);
                    }
            }
    }

    void UndirectedFirstTreeCalculator::SetDirections(boost::unordered_set<size_t> vertices_nums,
                                                      EvolutionaryTree& tree,
                                                      boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges) {
        auto edge_constructor = GetEdgeConstructor();
        auto undirected_graph = tree.GetUndirectedGraph();
        boost::unordered_set<size_t> undirected_graph_vertices;
        for (auto p : undirected_graph) {
            undirected_graph_vertices.insert(p.first);
        }

        for (auto clone_num : vertices_nums) {
            if (undirected_graph_vertices.find(clone_num) == undirected_graph_vertices.end()) {
                if (tree.GetUndirectedCompopentRoot(ds_on_undirected_edges.find_set(clone_num)) != size_t(-1)) {
                    const EvolutionaryEdge& edge = tree.GetUndirectedComponentParentEdge(clone_num);
                    tree.AddDirected(clone_num, edge, model_);
                };
                continue;
            }
            auto vertex_it = undirected_graph.find(clone_num);
            auto vertex = *vertex_it;
            std::vector<std::pair<size_t, size_t>> edge_vector;
            size_t root = tree.GetUndirectedCompopentRoot(ds_on_undirected_edges.find_set(vertex.first));
            if (root != size_t(-1)) {
                const EvolutionaryEdge& edge = tree.GetUndirectedComponentParentEdge(ds_on_undirected_edges.find_set(vertex.first));
                tree.AddDirected(edge.dst_clone_num, edge, model_);
                //tree.PrepareSubtree(edge_vector, root);

                //INFO("ds root is: " << ds_on_undirected_edges.find_set(vertex.first) << ", root is: " << root <<
                //" root's ds_root is:" << ds_on_undirected_edges.find_set(root));
                tree.PrepareSubtreeEdmonds(edge_vector, root, model_, clone_set_, edge_constructor);
            }
            else {
                //tree.PrepareSubtree(edge_vector, vertex.first);
                tree.PrepareSubtreeEdmonds(edge_vector, vertex.first, model_, clone_set_, edge_constructor);
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

    EvolutionaryEdgeConstructor* UndirectedFirstTreeCalculator::GetEdgeConstructor() {
        return new VJEvolutionaryEdgeConstructor(config_.edge_construction_params);
    }
}