#include <verify.hpp>
#include "vj_alignment_info.hpp"

namespace vj_finder {
    void VJAlignmentInfo::Update(VJAlignmentInfo vj_alignment) {
        VERIFY_MSG(false, "IMPLEMENT ME");
    }

    void VJAlignmentInfo::UpdateFilteredReads(const core::Read &read) {
        VERIFY_MSG(false, "IMPLEMENT ME");
    }

    void VJAlignmentInfo::UpdateHits(VJHits vj_hits) {
        VERIFY_MSG(false, "IMPLEMENT ME");
    }

    void VJAlignmentOutput::OutputAlignmentInfo() const {
        VERIFY_MSG(false, "IMPLEMENT ME");
    }

    void VJAlignmentOutput::OutputCleanedReads() const {
        VERIFY_MSG(false, "IMPLEMENT ME");
    }

    void VJAlignmentOutput::OutputFilteredReads() const {
        VERIFY_MSG(false, "IMPLEMENT ME");
    }
}