#pragma once

#include <boost/unordered_map.hpp>
#include "evolutionary_edge.hpp"
#include "model_utils/shm_model.hpp"

namespace antevolo {
    class EvolutionaryTree {
        boost::unordered_map<size_t, EvolutionaryEdge> edges_;
        boost::unordered_map<size_t, std::set<size_t>> undirected_graph_;
        boost::unordered_map<size_t, EvolutionaryEdge> undirected_components_edges_;

    public:
        void AddDirected(size_t clone_num, EvolutionaryEdge edge, ShmModel& model);

        void SetUndirectedComponentParentEdge(size_t root_num, EvolutionaryEdge edge, ShmModel& model);

        void AddUndirected(size_t clone_num, EvolutionaryEdge edge);

        void AddUndirectedPair(size_t src_num, size_t dst_num);

        void PrepareSubtree(std::vector<std::pair<size_t, size_t>>& edge_vector, size_t root_num);

        boost::unordered_map<size_t, std::set<size_t>>& GetUndirectedGraph()
        { return undirected_graph_; };

        bool Contains(size_t clone_num) {
            return (edges_.find(clone_num) != edges_.end());
        }

        const EvolutionaryEdge& GetParentEdge(size_t clone_num) {
            return edges_[clone_num];
        }

        size_t GetUndirectedCompopentRoot(size_t root_num) {
            if (undirected_components_edges_.find(root_num) != edges_.end()){
                return undirected_components_edges_[root_num].dst_clone_num;
            }
            else {
                return size_t(-1);
            }
        }

        const EvolutionaryEdge& GetUndirectedComponentParentEdge(size_t root_num) {
            return undirected_components_edges_[root_num];
        }

        void WriteInFile(std::string output_fname);

        size_t NumEdges() const { return edges_.size(); }

        size_t UndirectedGraphSize() const {
            size_t res = 0;
            for (auto it = undirected_graph_.begin(); it != undirected_graph_.end(); it++) {
                res += it->second.size();
            }
            return res;
        }
    };
}