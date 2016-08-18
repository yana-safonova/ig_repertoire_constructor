#include "tree_based_shm_annotator.hpp"

namespace antevolo {
    EvolutionaryAnnotatedSHM TreeBasedSHMAnnotator::GetAnnotation(size_t clone_id, annotation_utils::SHM shm) {
        VERIFY_MSG(false, "Implement me!");
        return EvolutionaryAnnotatedSHM(shm, false, false);
    }
}