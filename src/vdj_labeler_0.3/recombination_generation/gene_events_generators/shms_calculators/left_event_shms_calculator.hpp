#include "shm_calculator.hpp"

namespace vdj_labeler {

class LeftEventSHMsCalculator: public SHMsCalculator {
    int ComputeNumberCleavedSHMs(alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                 size_t left_cleavage_length);

    int ComputeNumberPalindromeSHMs(alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                    size_t left_palindrome_length);

public:
    int ComputeNumberSHMs(alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                          int left_cleavage_length,
                          int right_cleavage_length);

    int ComputeNumberSHMsForLeftEvent(alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                      int left_cleavage_length);

    int ComputeNumberSHMsForRightEvent(alignment_utils::ImmuneGeneReadAlignmentPtr, int) { return 0; }
};

} // End namespace vdj_labeler