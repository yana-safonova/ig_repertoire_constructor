#pragma once

#include "alignment_quality_checker.hpp"

namespace vdj_labeler {

class ThresholdAlignmentQualityChecker: public AlignmentQualityChecker {
    double normalized_score_threshold_;

public:
    ThresholdAlignmentQualityChecker(const double normalized_score_threshold) :
        normalized_score_threshold_(normalized_score_threshold) { }

    bool AlignmentIsGood(const alignment_utils::ImmuneGeneReadAlignment &ig_gene_alignment) const override {
        return ig_gene_alignment.NormalizedScore() >= normalized_score_threshold_;
    }
};

} // End namespace vdj_labeler