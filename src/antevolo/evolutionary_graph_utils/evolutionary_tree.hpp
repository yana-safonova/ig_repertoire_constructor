#pragma once

#include <boost/unordered_map.hpp>
#include "evolutionary_edge.hpp"

namespace antevolo {
    class EvolutionaryTree {
        boost::unordered_map<size_t, EvolutionaryEdge> edges_;
        boost::unordered_map<size_t, std::set<size_t>> undirected_graph_;

    public:
        void AddDirected(size_t clone_num, EvolutionaryEdge edge) {
            if(edge.IsDirected()) {
                if (!Contains(clone_num) || get_parent_edge(clone_num).num_added_v_shms > edge.num_added_v_shms) {
                    //if clone_set_[*it2] is root or if the new edge is shorter
                    edges_[clone_num] = edge;
                }
            }
        }
        void AddUndirected(size_t clone_num, EvolutionaryEdge edge) {
            if (edge.IsUndirected()) {
                edges_[clone_num] = edge;
            }
        }
        void AddUndirectedPair(size_t src_num, size_t dst_num) {
            if (undirected_graph_[src_num].find(dst_num) == undirected_graph_[src_num].end() &&
                    undirected_graph_[dst_num].find(src_num) == undirected_graph_[dst_num].end()) {
                undirected_graph_[src_num].insert(dst_num);
                undirected_graph_[dst_num].insert(src_num);
            }
        }

        void Prepare_subtree(std::vector<std::pair<size_t, size_t>>& edge_vector, size_t root_num) {
            if (undirected_graph_.find(root_num) != undirected_graph_.end()) {
                for (size_t u : undirected_graph_[root_num]) {
                    undirected_graph_[root_num].erase(u);
                    undirected_graph_[u].erase(root_num);
                    edge_vector.push_back(std::make_pair(root_num, u));
                    Prepare_subtree(edge_vector, u);
                }
            }
        }

        boost::unordered_map<size_t, std::set<size_t>>& Get_undirected_graph()
        { return undirected_graph_; };

        bool Contains(size_t clone_num) {
            return (edges_.find(clone_num) != edges_.end());
        }

        const EvolutionaryEdge& get_parent_edge(size_t clone_num) {
            return edges_[clone_num];
        }

        void WriteInFile(std::string output_fname);

        size_t NumEdges() const { return edges_.size(); }
    };
}