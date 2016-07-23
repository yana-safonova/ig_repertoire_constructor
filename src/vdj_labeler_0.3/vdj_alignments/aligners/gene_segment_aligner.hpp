#pragma once

#include "alignment_utils/alignment_positions.hpp"
#include "alignment_utils/pairwise_alignment.hpp"

namespace vdj_labeler {

class GeneSegmentAligner {
public:
    GeneSegmentAligner() { }
    GeneSegmentAligner(const GeneSegmentAligner &) = delete;
    GeneSegmentAligner& operator=(const GeneSegmentAligner &) = delete;
    GeneSegmentAligner(GeneSegmentAligner &&) = delete;
    GeneSegmentAligner& operator=(GeneSegmentAligner&&) = delete;

    virtual alignment_utils::ImmuneGeneReadAlignment ComputeAlignment(
        const alignment_utils::ImmuneGeneAlignmentPositions &alignment_positions) const = 0;
    virtual ~GeneSegmentAligner() { }
};

} // End namespace vdj_labeler