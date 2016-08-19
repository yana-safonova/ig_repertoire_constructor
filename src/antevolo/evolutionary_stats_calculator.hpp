#pragma once

#include "evolutionary_graph_utils/evolutionary_tree.hpp"
#include "evolutionary_tree_annotation/annotated_evolutionary_tree.hpp"

namespace antevolo {
    class EvolutionaryStatsCalculator {
        const annotation_utils::CDRAnnotatedCloneSet &clone_set_;

        AnnotatedEvolutionaryTree AddSHMsFromRoot(AnnotatedEvolutionaryTree annotated_tree);

        AnnotatedEvolutionaryTree AddSHMsFromEdges(AnnotatedEvolutionaryTree annotated_tree);

    public:
        EvolutionaryStatsCalculator(const annotation_utils::CDRAnnotatedCloneSet &clone_set) :
                clone_set_(clone_set) { }

        AnnotatedEvolutionaryTree ComputeTreeStats(const EvolutionaryTree &tree);
    };
}