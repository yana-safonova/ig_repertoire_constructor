#include "aa_shm_map.hpp"

namespace antevolo {
    bool operator<(const AaSHM &left, const AaSHM & right) {
        if(left.region != right.region)
            return left.region < right.region;
        if(left.consensus_pos != right.consensus_pos)
            return left.consensus_pos < right.consensus_pos;
        return left.dst_aa < right.dst_aa;
    }

    /*
    void AaSHMMap::AddSHM(TreeSHM tree_shm, size_t src_id, size_t dst_id) {
        VERIFY(false);
    }

    void AaSHMMap::DeleteEdge(size_t src_id, size_t dst_id) {
        VERIFY(false);
    }
    */
}