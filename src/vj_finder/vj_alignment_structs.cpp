#include "vj_alignment_structs.hpp"

#include <iostream>

namespace vj_finder {
    /*
    void VJHit::CheckConsistencyFatal() {
        VERIFY_MSG(v_hit_.Chain() == j_hit_.Chain(), "Chain of V gene hit (" << v_hit_.Chain() <<
                " ) does not match with chain of J gene hit (" << j_hit_.Chain() << ")");
    }

    void VJHit::AddLeftShift(int left_shift) {
        v_hit_.AddShift(left_shift);
        j_hit_.AddShift(left_shift);
    }

    void VJHit::AddRightShift(int right_shift) {
        j_hit_.AddShift(right_shift);
    }
     */

    std::vector<ImmuneGeneHitPtr> VJHits::VPtrHits() const {
        std::vector<ImmuneGeneHitPtr> v_ptr_hits;
        for (auto& hit : v_hits_) {
            v_ptr_hits.push_back(std::make_shared<VGeneHit>(hit));
        }
        return v_ptr_hits;
    }

    std::vector<ImmuneGeneHitPtr> VJHits::JPtrHits() const {
        std::vector<ImmuneGeneHitPtr> j_ptr_hits;
        for (auto& hit : j_hits_) {
            j_ptr_hits.push_back(std::make_shared<JGeneHit>(hit));
        }
        return j_ptr_hits;
    }
}