#include <clonally_related_candidates_calculators/edmonds_tarjan_DMST_calculator.hpp>
#include "evolutionary_tree.hpp"
#include "model_utils/shm_model.hpp"

namespace antevolo {
    void EvolutionaryTree::AddDirected(size_t clone_num, EvolutionaryEdge edge, ShmModel& model) {
        if(edge.IsDirected()) {
            if (!Contains(clone_num)) {
                edges_[clone_num] = edge;
                return;
            }
            const EvolutionaryEdge& parent_edge =  edges_[clone_num];
            if (static_cast<double>(parent_edge.num_added_shms) + 1.0 - model.CDR3TransitionProb(parent_edge) >
                static_cast<double>(edge.num_added_shms) + 1.0 - model.CDR3TransitionProb(edge)) {
                //if clone_set_[*it2] is root or if the new edge is shorter
                edges_[clone_num] = edge;
            }
        }
    }

    void EvolutionaryTree::SetUndirectedComponentParentEdge(size_t root_num, EvolutionaryEdge edge, ShmModel& model) {
        if(edge.IsDirected()) {
            if (undirected_components_edges_.find(root_num) == undirected_components_edges_.end()) {
                undirected_components_edges_[root_num] = edge;
                return;
            }
            const EvolutionaryEdge& parent_edge =  undirected_components_edges_[root_num];
            if (static_cast<double>(parent_edge.num_added_shms) + 1.0 - model.CDR3TransitionProb(parent_edge) >
                static_cast<double>(edge.num_added_shms) + 1.0 - model.CDR3TransitionProb(edge)) {
                //if clone_set_[*it2] is root or if the new edge is shorter
                undirected_components_edges_[root_num] = edge;
            }
        }
    }

    void EvolutionaryTree::AddUndirected(size_t clone_num, EvolutionaryEdge edge) {
        if (edge.IsUndirected()) {
            edges_[clone_num] = edge;
        }
    }

    void EvolutionaryTree::AddUndirectedPair(size_t src_num, size_t dst_num) {
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

    void EvolutionaryTree::PrepareSubtree(std::vector<std::pair<size_t, size_t>>& edge_vector, size_t root_num) {
        if (undirected_graph_.find(root_num) != undirected_graph_.end()) {
            for (size_t u : undirected_graph_[root_num]) {
                undirected_graph_[root_num].erase(u);
                undirected_graph_[u].erase(root_num);
                edge_vector.push_back(std::make_pair(root_num, u));
                PrepareSubtree(edge_vector, u);
            }
        }
    }

    void EvolutionaryTree::PrepareSubtreeVertices(
            boost::unordered_set<size_t>& vertices_set,
            size_t root_num) {
        //we assume that component size is > 1
        vertices_set.insert(root_num);
        for (size_t u : undirected_graph_[root_num]) {
            undirected_graph_[root_num].erase(u);
            undirected_graph_[u].erase(root_num);
            PrepareSubtreeVertices(vertices_set, u);
        }
    }

    void EvolutionaryTree::PrepareSubtreeEdmonds(std::vector<std::pair<size_t, size_t>>& edge_vector,
                                                 size_t root_vertex,
                                                 ShmModel& model,
                                                 const annotation_utils::CDRAnnotatedCloneSet& clone_set,
                                                 EvolutionaryEdgeConstructor* edge_constructor) {
        //INFO("Preparing graph, root is " << root_vertex);
        boost::unordered_set<size_t> vertices_set;
        PrepareSubtreeVertices(vertices_set, root_vertex);
        for (size_t v : vertices_set) {
            undirected_graph_.erase(v);
        }
        size_t n = vertices_set.size();
        boost::unordered_map<size_t, size_t> vertex_to_index;

        /*
        for (auto v : vertices_set) {
            std::cout << v << " ";
        } std::cout << std::endl;
        */
        std::vector<size_t> index_to_vertex(n);
        size_t index = 0;
        for (size_t vertex : vertices_set) {
            //std::cout << clone_to_string(clone_set[vertex]) << std::endl << std::endl;
            vertex_to_index[vertex] = index;
            index_to_vertex[index] = vertex;
            ++index;
        }
        size_t root_index = vertex_to_index[root_vertex];

        typedef EdmondsTarjanDMSTCalculator::WeightedEdge WeightedEdge;
        std::vector<WeightedEdge> edges;
        for (size_t v = 0; v < n; ++v) {
            for (size_t u = 0; u < n; ++u) {
                if (v == u) {
                    continue;
                }
                //INFO(index_to_vertex[u] << "->" << index_to_vertex[v]);
                double dist = model.CDR3TransitionProb(
                        edge_constructor->ConstructEdge(clone_set[index_to_vertex[v]],
                                                        clone_set[index_to_vertex[u]],
                                                        index_to_vertex[v],
                                                        index_to_vertex[u]));
                if (dist != 0) {
                    edges.push_back(WeightedEdge(v, u, log(dist)));
                    //std::cout << v << " " << u << " " << log(dist)  << std::endl;
                }

            }
        }
        //std::cout << "root " << root_index  << std::endl;
        //INFO("Graph prepared");
        EdmondsTarjanDMSTCalculator e_calc(n, edges, root_index);
        e_calc.EmpondsTarjan();
        std::vector<WeightedEdge> edges_to_add = e_calc.GetParentEdges();
        for (auto& we : edges_to_add) {
            size_t src_clone_num = index_to_vertex[we.src_];
            size_t dst_clone_num = index_to_vertex[we.dst_];
            edge_vector.push_back(std::make_pair(src_clone_num, dst_clone_num));
            //std::cout << we.src_ << " " << we.dst_ << " (" << we.weight_ <<  ")"  << std::endl;
        }
        //INFO("Minumal arborescence found");
    }

    void EvolutionaryTree::WriteEdge(const EvolutionaryEdge& edge, std::ofstream& out) { //no endl
        out << edge.src_clone->Read().id << "\t" << edge.dst_clone->Read().id << "\t"
        << edge.src_clone->Read().name << "\t" << edge.dst_clone->Read().name << "\t"
        << edge.src_clone->VSHMs().size() + edge.src_clone->JSHMs().size() << "\t"
        << edge.dst_clone->VSHMs().size() + edge.dst_clone->JSHMs().size()<< "\t"
        << edge.edge_type << "\t" <<  edge.num_intersected_shms << "\t" << edge.num_added_shms
        << "\t" << edge.cdr3_distance << "\t" << edge.weight << "\t";
    }

    void EvolutionaryTree::WriteInFile(std::string output_fname) {
        std::ofstream out(output_fname);
        out << "Src_id\tDst_id\tSrc_clone\tDst_clone\tNum_Src_SHMs\tNum_Dst_SHMs\tEdge_type\tNum_shared_SHMs\tNum_added_SHMs\tCDR3_dist\tWeight\n";
        for(auto it = edges_.begin(); it != edges_.end(); it++)
            WriteEdge(it->second, out);
        out << std::endl;
        out.close();
    }

    void EvolutionaryTree::WriteInFileWithCDR3s(std::string output_fname) {
        std::ofstream out(output_fname);
        out <<
        "Src_id\tDst_id\tSrc_clone\tDst_clone\tNum_Src_SHMs\tNum_Dst_SHMs\tEdge_type\tNum_shared_SHMs\tNum_added_SHMs\tCDR3_dist\tWeight\tSrc_CDR3\tDst_CDR3\n";
        for(auto it = edges_.begin(); it != edges_.end(); it++) {
            WriteEdge(it->second, out);
            std::string src_CDR3_string;
            std::string dst_CDR3_string;
            auto const& src_CDR3 = it->second.src_clone->CDR3();
            auto const& dst_CDR3 = it->second.dst_clone->CDR3();
            size_t src_CDR3_length = seqan::length(src_CDR3);
            size_t dst_CDR3_length = seqan::length(dst_CDR3);
            for (size_t pos = 0; pos < src_CDR3_length; ++pos) {
                src_CDR3_string.push_back(src_CDR3[pos]);
            }
            for (size_t pos = 0; pos < dst_CDR3_length; ++pos) {
                dst_CDR3_string.push_back(dst_CDR3[pos]);
            }
            out << src_CDR3_string << "\t" << dst_CDR3_string << std::endl;
        }
        out.close();
    }
}