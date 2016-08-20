#pragma once

#include "evolutionary_tree.hpp"

namespace antevolo {
    class EvolutionaryTreeSplitter {
    public:
        virtual std::vector<EvolutionaryTree> Split(const EvolutionaryTree &tree) = 0;

        virtual ~EvolutionaryTreeSplitter() { }
    };

    class ConnectedTreeSplitter : public EvolutionaryTreeSplitter {
    public:
        std::vector<EvolutionaryTree> Split(const EvolutionaryTree &tree);
    };
}