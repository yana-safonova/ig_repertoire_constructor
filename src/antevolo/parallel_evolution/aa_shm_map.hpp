#pragma once

#include "../shm_counting/tree_based_shm.hpp"

namespace antevolo {
    struct AaSHM {
        annotation_utils::StructuralRegion region;
        size_t consensus_pos;
        char dst_aa;

        AaSHM(annotation_utils::StructuralRegion region, size_t consensus_pos,
              char dst_aa) : region(region),
                             consensus_pos(consensus_pos),
                             dst_aa(dst_aa) { }
    };

    bool operator<(const AaSHM &left, const AaSHM & right);

    class AaSHMMap {
        typedef std::pair<size_t, size_t> ClonalEdge;
        std::map<AaSHM, size_t> aa_mult_map_;
        std::map<AaSHM, std::vector<ClonalEdge> > aa_dst_id_map_;
        std::map<ClonalEdge, std::vector<AaSHM> > edge_shms_map_;

    public:
        AaSHMMap() { }

        /*
        void AddSHM(TreeSHM tree_shm, size_t src_id, size_t dst_id);

        void DeleteEdge(size_t src_id, size_t dst_id);
        */
    };
}