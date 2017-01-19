#pragma once

#include "evolutionary_graph_utils/evolutionary_tree.hpp"
#include "evolutionary_tree_annotation/annotated_tree_storage.hpp"
#include "evolutionary_tree_storage.hpp"

namespace antevolo {
    class EvolutionaryStatsCalculator {

        AnnotatedEvolutionaryTree AddSHMsFromRoot(AnnotatedEvolutionaryTree annotated_tree) const;

        AnnotatedEvolutionaryTree AddSHMsFromEdges(AnnotatedEvolutionaryTree annotated_tree) const;

    public:

        AnnotatedEvolutionaryTree ComputeStatsForTree(const EvolutionaryTree &tree) const;

        AnnotatedTreeStorage ComputeStatsForStorage(const EvolutionaryTreeStorage& storage) const;
    };
}