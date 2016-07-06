#pragma once

#include "alignment_quality_checker.hpp"

namespace vdj_labeler {

class ThresholdAlignmentQualityChecker: public AlignmentQualityChecker {
    double normalized_score_threshold_;

public:
    ThresholdAlignmentQualityChecker(double normalized_score_threshold) :
        normalized_score_threshold_(normalized_score_threshold) { }

    bool AlignmentIsGood(alignment_utils::ImmuneGeneReadAlignmentPtr ig_gene_alignment) {
        return ig_gene_alignment->NormalizedScore() >= normalized_score_threshold_;
    }
};

} // End namespace vdj_labeler