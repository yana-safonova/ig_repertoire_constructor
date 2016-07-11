#pragma once

#include <boost/unordered_map.hpp>
#include "evolutionary_edge.hpp"

namespace antevolo {
    class EvolutionaryTree {
        boost::unordered_map<size_t, EvolutionaryEdge> edges_;

    public:
        void Add(size_t clone_num, EvolutionaryEdge edge) {
            if(edge.IsDirected()) {
                if (!Contains(clone_num) || Get_parent_edge(clone_num).num_added_v_shms > edge.num_added_v_shms) {
                    //if clone_set_[clone-num] is root or if the new edge is shorter
                    edges_[clone_num] = edge;
                }
            }
        }
        void AddUndirected(size_t clone_num, EvolutionaryEdge edge) {
            if(edge.IsUndirected()) {
                if (!Contains(clone_num)) {
                    //if clone_set_[clone_num] is root
                    edges_[clone_num] = edge;
                }
            }
        }

        bool Contains(size_t clone_num) {
            return (edges_.find(clone_num) != edges_.end());
        }

        const EvolutionaryEdge& Get_parent_edge(size_t clone_num) {
            return edges_[clone_num];
        }

        void WriteInFile(std::string output_fname);

        size_t NumEdges() const { return edges_.size(); }
    };
}