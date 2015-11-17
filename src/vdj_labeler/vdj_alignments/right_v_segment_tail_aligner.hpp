#pragma once

#include "gene_segment_aligner.hpp"

class RightVSegmentTailAligner : public GeneSegmentAligner {
public:
    IgGeneAlignmentPtr ComputeAlignment(IgGeneAlignmentPositions alignment_positions);
};
