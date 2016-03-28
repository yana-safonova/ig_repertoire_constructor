#include "pairwise_block_alignment.hpp"

namespace algorithms {
    PairwiseBlockAlignment PairwiseBlockAlignment::path2Alignment(AlignmentPath &path, const seqan::Dna5String &read,
                                                                  const seqan::Dna5String &query, size_t needle_index,
                                                                  int score) {
        int coverage_length = path.kplus_length();

        int left_shift = path.left_shift();
        int right_shift = path.right_shift();

        int start = left_shift;
        int finish = right_shift + int(seqan::length(query));

        int over_start = std::max(0, start);
        int over_finish = std::min(right_shift + seqan::length(query), seqan::length(read));
        int read_overlap_length = over_finish - over_start; // read overlap
        int needle_overlap_length = read_overlap_length + left_shift - right_shift;

        PairwiseBlockAlignment align;
        align.kp_coverage = coverage_length;
        align.int_score = score;
        // align.score = static_cast<double>(coverage_length) / static_cast<double>(length(query));
        align.score = static_cast<double>(score) / static_cast<double>(seqan::length(query));
        align.path = std::move(path);
        align.start = start;
        align.finish = finish;
        align.needle_index = needle_index;
        align.overlap_length = needle_overlap_length;
        align.needle_length = seqan::length(query);

        return align;

    }
}