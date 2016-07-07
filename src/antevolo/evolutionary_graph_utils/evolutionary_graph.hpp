#pragma once

#include "evolutionary_edge.hpp"
#include "evolutionary_tree.hpp"

namespace antevolo {
    struct AnnotatedCloneHasher {
        size_t operator()(const annotation_utils::AnnotatedClone& clone) const {
            return std::hash<std::string>()(clone.Read().name);
        }
    };

    struct EvolutionaryEdgeHasher {
        size_t operator()(const EvolutionaryEdge& edge) const {
            return AnnotatedCloneHasher().operator()(*edge.src_clone) *
                   AnnotatedCloneHasher().operator()(*edge.dst_clone);
        }
    };

    class EvolutionaryGraph {
        std::vector<EvolutionaryEdge> edges_;
        //std::unordered_map<EvolutionaryEdge, size_t, EvolutionaryEdgeHasher> edge_index_map_;
        std::unordered_map<annotation_utils::AnnotatedClone, size_t, AnnotatedCloneHasher> vertex_index_map_;

        std::map<size_t, std::vector<size_t>> vertex_edges_map_;
        std::map<size_t, std::vector<size_t>> vertex_neighbours_map_;

        std::unordered_map<EvolutionaryEdgeType, size_t, std::hash<int>> edge_type_counts_;

        // if it is not presented in map, add new record
        size_t GetVertexIndex(const annotation_utils::AnnotatedClone& clone);

        void UpdateVertexMap(size_t src_index, size_t dst_index);

        void UpdateEdgeMap(size_t vertex_index, size_t edge_index);

    public:
        void AddEdge(EvolutionaryEdge edge);

        size_t NumEdges() const { return edges_.size(); }

        size_t NumEdgesOfType(EvolutionaryEdgeType edge_type) const;

        EvolutionaryTree GetEvolutionaryTree();

        typedef std::vector<EvolutionaryEdge>::const_iterator EdgeIterator;

        EdgeIterator cbegin() const { return edges_.cbegin(); }

        EdgeIterator cend() const { return edges_.cend(); }
    };
}