#pragma once

#include "../alignment_structs.hpp"

class GeneSegmentAligner {
public:
    virtual IgGeneAlignmentPtr ComputeAlignment(IgGeneAlignmentPositions alignment_positions) = 0;
    virtual ~GeneSegmentAligner() { }
};