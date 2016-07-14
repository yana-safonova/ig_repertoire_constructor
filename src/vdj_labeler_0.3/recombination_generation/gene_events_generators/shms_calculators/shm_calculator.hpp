#pragma once

#include "alignment_utils/pairwise_alignment.hpp"

namespace vdj_labeler {

class SHMsCalculator {
public:
    virtual int ComputeNumberSHMs(const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                  const int left_cleavage_length,
                                  const int right_cleavage_length) const = 0;

    virtual int ComputeNumberSHMsForLeftEvent(const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                              const int left_cleavage_length) const = 0;

    virtual int ComputeNumberSHMsForRightEvent(const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                               const int right_cleavage_length) const = 0;

    virtual ~SHMsCalculator() { }
};

} // End namespace vdj_labeler