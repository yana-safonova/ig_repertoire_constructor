#include "tree_based_shm_convertor.hpp"

namespace antevolo {
    size_t TreeSHMComparator::GetTreeSHMPosition(const annotation_utils::AnnotatedClone &related_clone, TreeSHM shm) {
        if(shm.region == annotation_utils::CDR3)
            return related_clone.CDR3Range().start_pos + shm.gene_pos;
        if(shm.region == annotation_utils::FR4) {
            return related_clone.JAlignment().QueryPositionBySubjectPosition(shm.gene_pos);
        }
        return related_clone.VAlignment().QueryPositionBySubjectPosition(shm.gene_pos);
    }
}