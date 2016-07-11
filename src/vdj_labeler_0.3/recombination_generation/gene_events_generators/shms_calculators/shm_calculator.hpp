#pragma once

#include "alignment_utils/pairwise_alignment.hpp"

namespace vdj_labeler {

class SHMsCalculator {
public:
    virtual int ComputeNumberSHMs(alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                  int left_cleavage_length,
                                  int right_cleavage_length) = 0;

    virtual int ComputeNumberSHMsForLeftEvent(alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                              int left_cleavage_length) = 0;

    virtual int ComputeNumberSHMsForRightEvent(alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                               int right_cleavage_length) = 0;

    virtual ~SHMsCalculator() { }
};

} // End namespace vdj_labeler