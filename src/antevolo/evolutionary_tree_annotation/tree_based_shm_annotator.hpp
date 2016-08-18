#pragma once

#include "evolutionary_graph_utils/evolutionary_tree.hpp"
#include "evolutionary_annotated_shm.hpp"

namespace antevolo {
    class TreeBasedSHMAnnotator {
        const annotation_utils::CDRAnnotatedCloneSet& clone_set_;
        const EvolutionaryTree &tree_;

    public:
        TreeBasedSHMAnnotator(const annotation_utils::CDRAnnotatedCloneSet& clone_set,
                              const EvolutionaryTree &tree) : clone_set_(clone_set),
                                                              tree_(tree) { }

        EvolutionaryAnnotatedSHM GetAnnotation(size_t clone_id, annotation_utils::SHM shm);
    };
}