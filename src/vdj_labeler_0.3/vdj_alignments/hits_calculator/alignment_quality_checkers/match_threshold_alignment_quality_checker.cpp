#include "match_threshold_alignment_estimator.hpp"

bool MatchThresholdAlignmentEstimator::AlignmentIsGood(IgGeneAlignmentPtr ig_gene_alignment) {
    int cnt_following_matches = 0;
    auto row0 = seqan::row(ig_gene_alignment->Alignment(), 0);
    auto row1 = seqan::row(ig_gene_alignment->Alignment(), 1);
    for (auto it_row0 = begin(row0), it_row1 = begin(row1); it_row0 != end(row0); ++it_row0, ++it_row1) {
        if (*it_row0 == *it_row1)
            cnt_following_matches++;
        else {
            cnt_following_matches = 0;
            continue;
        }
        if (cnt_following_matches == this->normalized_score_threshold_)
            return true;
    }
    return false;
}
