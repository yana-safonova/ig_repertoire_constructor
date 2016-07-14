#pragma once

#include "shm_calculator.hpp"

namespace vdj_labeler {

class RightEventSHMsCalculator : public SHMsCalculator {
    int ComputeNumberCleavedSHMs(const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                 const size_t cleavage_length) const;

    int ComputeNumberPalindromeSHMs(const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                    const size_t palindrome_length) const;

public:
    int ComputeNumberSHMs(const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                          const int left_cleavage_length,
                          const int right_cleavage_length) const;

    int ComputeNumberSHMsForLeftEvent(const alignment_utils::ImmuneGeneReadAlignmentPtr, const int) const { return 0; }

    int ComputeNumberSHMsForRightEvent(const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                       const int right_cleavage_length) const;
};

} // End namespace vdj_labeler