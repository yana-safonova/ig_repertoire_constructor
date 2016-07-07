#pragma once

#include "evolutionary_edge.hpp"

namespace antevolo {
    class EvolutionaryTree {
        std::vector<EvolutionaryEdge> edges_;

    public:
        void Add(EvolutionaryEdge edge) {
            edges_.push_back(edge);
        }

        void WriteInFile(std::string output_fname);

        size_t NumEdges() const { return edges_.size(); }
    };
}