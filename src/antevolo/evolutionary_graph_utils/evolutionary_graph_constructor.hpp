#pragma once

#include "evolutionary_edge_constructor.hpp"
#include "evolutionary_graph.hpp"
#include "../clonally_related_candidates_calculators/clonally_related_candidates.hpp"

namespace antevolo {
    class EvolutionaryGraphConstructor {
        const AntEvoloConfig::AlgorithmParams &config_;
        const annotation_utils::CDRAnnotatedCloneSet &clone_set_;
        const ClonallyRelatedCandidates candidates_;

        EvolutionaryEdgeConstructor* GetEdgeConstructor();

    public:
        EvolutionaryGraphConstructor(const AntEvoloConfig::AlgorithmParams &config,
                                     const annotation_utils::CDRAnnotatedCloneSet &clone_set,
                                     const ClonallyRelatedCandidates candidates) : config_(config),
                                                                                   clone_set_(clone_set),
                                                                                   candidates_(candidates) { }

       // EvolutionaryGraph ConstructGraph();
    };
}