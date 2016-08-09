//
// Created by Andrew Bzikadze on 8/5/16.
//


#include "logger/logger.hpp"

#include "single_immune_gene_segment_hits.hpp"
#undef NDEBUG
#include <cassert>

using namespace vdj_labeler;

void SingleImmuneGeneSegmentHits::AddHit(alignment_utils::ImmuneGeneReadAlignment hit) {
    if (hit.SubjectPtr() != nullptr) {
        VERIFY(hit.Subject().GeneType().Segment() == segment_type_);
    }
    hits_.emplace_back(std::move(hit));
}

const alignment_utils::ImmuneGeneReadAlignment& SingleImmuneGeneSegmentHits::operator[](const size_t &index) const {
    VERIFY(index < size());
    return hits_[index];
}
