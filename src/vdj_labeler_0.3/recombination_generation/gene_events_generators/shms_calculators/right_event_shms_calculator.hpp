#pragma once

#include "shm_calculator.hpp"

namespace vdj_labeler {

class RightEventSHMsCalculator : public SHMsCalculator {
    int ComputeNumberCleavedSHMs(alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                 size_t cleavage_length);

    int ComputeNumberPalindromeSHMs(alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                    size_t palindrome_length);

public:
    int ComputeNumberSHMs(alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                          int left_cleavage_length,
                          int right_cleavage_length);

    int ComputeNumberSHMsForLeftEvent(alignment_utils::ImmuneGeneReadAlignmentPtr, int) { return 0; }

    int ComputeNumberSHMsForRightEvent(alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                       int right_cleavage_length);
};

} // End namespace vdj_labeler