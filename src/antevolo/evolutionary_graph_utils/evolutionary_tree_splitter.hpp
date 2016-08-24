#pragma once

#include "evolutionary_tree.hpp"

namespace antevolo {
    class EvolutionaryTreeSplitter {

    public:
        virtual std::vector<EvolutionaryTree> Split(const EvolutionaryTree &tree) = 0;

        virtual ~EvolutionaryTreeSplitter() { }
    };

    class ConnectedTreeSplitter : public EvolutionaryTreeSplitter {
        EvolutionaryTree GetTreeByRoot(const EvolutionaryTree &tree, size_t root_id);

    public:
        std::vector<EvolutionaryTree> Split(const EvolutionaryTree &tree);
    };
}