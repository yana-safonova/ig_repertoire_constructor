#include "vj_alignment_structs.hpp"

namespace vj_finder {
    void VJHit::CheckConsistency() {
        VERIFY_MSG(v_hit_.Chain() == j_hit_.Chain(), "Chain of V gene hit (" << v_hit_.Chain() <<
                " ) does not match with chain of J gene hit (" << j_hit_.Chain() << ")");
    }
}