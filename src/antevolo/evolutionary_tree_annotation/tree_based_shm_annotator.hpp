#pragma once

#include "evolutionary_graph_utils/evolutionary_tree.hpp"
#include "evolutionary_annotated_shm.hpp"

namespace antevolo {
    class TreeBasedSHMAnnotator {
        const EvolutionaryTree *tree_ptr_;

        EvolutionaryAnnotatedSHM GetAnnotationForRoot(size_t root_id, annotation_utils::SHM shm);

        EvolutionaryAnnotatedSHM GetAnnotationForNodeWithParent(size_t clone_id, annotation_utils::SHM shm);

        bool SHMIsSynonymousWrtToParent(size_t clone_id, annotation_utils::SHM shm);

    public:
        TreeBasedSHMAnnotator(const EvolutionaryTree &tree) :
                tree_ptr_(&tree) { }

        EvolutionaryAnnotatedSHM GetAnnotation(size_t clone_id, annotation_utils::SHM shm);
    };
}