#pragma once

#include "gene_segment_aligner.hpp"

class LeftJTailAligner : public GeneSegmentAligner {
    size_t left_shift_;

public:
    LeftJTailAligner(size_t left_shift = 0) :
        left_shift_(left_shift) { }

    IgGeneAlignmentPtr ComputeAlignment(IgGeneAlignmentPositions alignment_positions);
};