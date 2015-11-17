#pragma once

#include "gene_segment_aligner.hpp"

class LeftJTailAligner : public GeneSegmentAligner {
public:
    IgGeneAlignmentPtr ComputeAlignment(IgGeneAlignmentPositions alignment_positions);
};