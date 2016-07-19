#pragma once

#include "alignment_quality_checker.hpp"

namespace vdj_labeler {

class MatchThresholdAlignmentQualityChecker: public AlignmentQualityChecker {
    int normalized_score_threshold_;

public:
    /*
     * The default value is 5 following IgBlast
     * link: https://nar.oxfordjournals.org/content/early/2013/05/11/nar.gkt382.full.pdf :
     * Identifying the V, D and J gene hits:
     * "The default word size is 5, which requires a minimum of five consecutive nucleotide
     * matches for a D gene to be found."
     */
    MatchThresholdAlignmentQualityChecker(int normalized_score_threshold = 5) :
        normalized_score_threshold_(normalized_score_threshold) { }

    bool AlignmentIsGood(alignment_utils::ImmuneGeneReadAlignmentPtr ig_gene_alignment) const override;
};

} // End namespace vdj_labeler
