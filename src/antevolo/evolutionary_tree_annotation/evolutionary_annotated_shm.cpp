#include "evolutionary_annotated_shm.hpp"

namespace antevolo {
    std::ostream& operator<<(std::ostream &out, const EvolutionaryAnnotatedSHM &shm) {
        out << shm.shm << ", synonymous: " << shm.synonymous << ", fixed: " << shm.fixed;
        return out;
    }

    bool CDR3SHMComparator::operator()(const CDR3SHM &shm1, const CDR3SHM &shm2) {
        if (shm1.rel_dst_pos != shm2.rel_dst_pos)
            return shm1.rel_dst_pos < shm2.rel_dst_pos;
        if (shm1.rel_src_pos != shm2.rel_src_pos)
            return shm1.rel_src_pos < shm2.rel_src_pos;
        if (shm1.dst_nucl != shm2.dst_nucl)
            return shm1.dst_nucl < shm2.dst_nucl;
        if (shm1.src_nucl != shm2.src_nucl)
            return shm1.src_nucl < shm2.src_nucl;
        if (shm1.dst_aa != shm2.dst_aa)
            return shm1.dst_aa < shm2.dst_aa;
        return shm1.src_aa < shm2.src_aa;
    }
}