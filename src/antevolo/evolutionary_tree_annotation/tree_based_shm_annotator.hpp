#pragma once

#include "evolutionary_graph_utils/evolutionary_tree.hpp"
#include "evolutionary_annotated_shm.hpp"

namespace antevolo {
    class TreeBasedSHMAnnotator {
        const annotation_utils::CDRAnnotatedCloneSet* clone_set_ptr_;
        const EvolutionaryTree *tree_ptr_;

        EvolutionaryAnnotatedSHM GetAnnotationForRoot(size_t root_id, annotation_utils::SHM shm);

        EvolutionaryAnnotatedSHM GetAnnotationForNodeWithParent(size_t clone_id, annotation_utils::SHM shm);

        bool SHMIsSynonymousWrtToParent(size_t clone_id, annotation_utils::SHM shm);

    public:
        TreeBasedSHMAnnotator(const annotation_utils::CDRAnnotatedCloneSet& clone_set,
                              const EvolutionaryTree &tree) : clone_set_ptr_(&clone_set),
                                                              tree_ptr_(&tree) { }

        EvolutionaryAnnotatedSHM GetAnnotation(size_t clone_id, annotation_utils::SHM shm);
    };
}