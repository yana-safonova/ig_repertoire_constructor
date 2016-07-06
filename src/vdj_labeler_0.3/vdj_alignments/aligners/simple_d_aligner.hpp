#pragma once

#include "gene_segment_aligner.hpp"

namespace vdj_labeler {

class SimpleDAligner {
public:
    alignment_utils::ImmuneGeneReadAlignmentPtr ComputeAlignment(
        alignment_utils::ImmuneGeneAlignmentPositions alignment_positions);
};

} // End namespace vdj_labeler
