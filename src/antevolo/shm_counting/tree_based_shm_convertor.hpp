#pragma once

#include "tree_shm_map.hpp"

namespace antevolo {
    class TreeSHMComparator {
    public:
        // get position of SHM on related clone
        // based on information about V/J alignment or CDR3 in clone with SHM
        static size_t GetTreeSHMPosition(const annotation_utils::AnnotatedClone &related_clone, TreeSHM shm);
    };
}