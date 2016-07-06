#pragma once

#include "alignment_quality_controller.hpp"

class ThresholdAlignmentEstimator : public AlignmentEstimator {
    double normalized_score_threshold_;

public:
    ThresholdAlignmentEstimator(double normalized_score_threshold) :
            normalized_score_threshold_(normalized_score_threshold) { }

    bool AlignmentIsGood(IgGeneAlignmentPtr ig_gene_alignment) {
        return ig_gene_alignment->NormalizedScore() >= normalized_score_threshold_;
    }
};