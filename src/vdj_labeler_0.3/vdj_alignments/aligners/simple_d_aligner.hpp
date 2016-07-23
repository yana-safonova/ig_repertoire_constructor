#pragma once

#include "gene_segment_aligner.hpp"

namespace vdj_labeler {

class SimpleDAligner : public GeneSegmentAligner {
public:
    virtual alignment_utils::ImmuneGeneReadAlignment ComputeAlignment(
        const alignment_utils::ImmuneGeneAlignmentPositions &alignment_positions) const override;
};

} // End namespace vdj_labeler
