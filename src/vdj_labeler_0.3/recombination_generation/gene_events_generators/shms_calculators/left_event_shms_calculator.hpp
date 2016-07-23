#pragma once

#include "shm_calculator.hpp"

namespace vdj_labeler {

class LeftEventSHMsCalculator: public SHMsCalculator {
    int ComputeNumberCleavedSHMs(const alignment_utils::ImmuneGeneReadAlignment& gene_alignment,
                                 const size_t left_cleavage_length) const;

    int ComputeNumberPalindromeSHMs(const alignment_utils::ImmuneGeneReadAlignment& gene_alignment,
                                    const size_t left_palindrome_length) const;

public:
    int ComputeNumberSHMs(const alignment_utils::ImmuneGeneReadAlignment& gene_alignment,
                          const int left_cleavage_length,
                          const int right_cleavage_length) const override;

    int ComputeNumberSHMsForLeftEvent(const alignment_utils::ImmuneGeneReadAlignment& gene_alignment,
                                      const int left_cleavage_length) const override;

    int ComputeNumberSHMsForRightEvent(const alignment_utils::ImmuneGeneReadAlignment&, const int) const override
    { return 0; }
};

} // End namespace vdj_labeler