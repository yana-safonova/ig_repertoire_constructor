#pragma once

#include "alignment_utils/alignment_positions.hpp"
#include "alignment_utils/pairwise_alignment.hpp"

namespace vdj_labeler {

class GeneSegmentAligner {
public:
    virtual alignment_utils::ImmuneGeneReadAlignmentPtr ComputeAlignment(
        alignment_utils::ImmuneGeneAlignmentPositions alignment_positions) = 0;
    virtual ~GeneSegmentAligner() { }
};

} // End namespace vdj_labeler