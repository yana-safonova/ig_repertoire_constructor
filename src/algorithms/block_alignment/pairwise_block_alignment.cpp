#include "pairwise_block_alignment.hpp"

namespace algorithms {
    PairwiseBlockAlignment::PairwiseBlockAlignment(AlignmentPath &path,
                                                   size_t subject_length,
                                                   size_t query_length,
                                                   int score) {
        this->path = std::move(path);
        this->subject_length = subject_length;
        this->query_length = query_length;

        kp_coverage = path.kplus_length();
        int_score = score;
        this->score = static_cast<double>(score) / static_cast<double>(subject_length);

        int left_shift = path.left_shift();
        int right_shift = path.right_shift();
        start = left_shift;
        finish = right_shift + int(subject_length);
        int over_start = std::max<int>(0, start);
        int over_finish = std::min<int>(right_shift + subject_length, static_cast<int>(query_length));
        int read_overlap_length = over_finish - over_start; // read overlap
        overlap_length = read_overlap_length + left_shift - right_shift;
    }
}