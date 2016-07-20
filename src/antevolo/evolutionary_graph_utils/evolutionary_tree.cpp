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
            if (parent_edge.num_added_v_shms + 1 - model.CDR3TransitionProb(parent_edge) >
                edge.num_added_v_shms + 1 - model.CDR3TransitionProb(edge)) {
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
            if (parent_edge.num_added_v_shms + 1 - model.CDR3TransitionProb(parent_edge) >
                edge.num_added_v_shms + 1 - model.CDR3TransitionProb(edge)) {
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

    void EvolutionaryTree::WriteInFile(std::string output_fname) {
        std::ofstream out(output_fname);
        out << "Src_num\tDst_num\tSrc_clone\tDst_clone\tNum_Src_V_SHMs\tNum_Dst_V_SHMs\tEdge_type\tNum_shared_SHMs\tNum_added_SHMs\tCDR3_dist\tWeight\n";
        for(auto it = edges_.begin(); it != edges_.end(); it++)
            out << it->second.src_clone_num << "\t" << it->second.dst_clone_num << "\t" << it->second.src_clone->Read().name << "\t" << it->second.dst_clone->Read().name << "\t" <<
            it->second.src_clone->VSHMs().size() << "\t" << it->second.dst_clone->VSHMs().size() << "\t" << it->second.edge_type << "\t" <<
            it->second.num_intersected_v_shms << "\t" << it->second.num_added_v_shms << "\t" << it->second.cdr3_distance << "\t" <<
            it->second.weight << std::endl;
        out.close();
    }
}