#pragma once

#include "verify.hpp"

#include "block_alignment_primitives.hpp"
#include "../hashes/subject_query_kmer_index.hpp"

#include <seqan/seq_io.h>
#include <seqan/align.h>

namespace algorithms {
    template<typename Tsequence1, typename Tsequence2>
    int find_simple_gap(const Tsequence1 &read, const Tsequence2 &gene) {
        // Positive if gaps in read (|read| < |gene|)
        // ggggggggg
        // rrrrrr
        //    rrrrrr
        using seqan::length;
        VERIFY(length(read) < length(gene));

        std::vector<int> cum_matches_forward(length(read) + 1), cum_matches_backward(length(read) + 1);

        for (size_t i = 1; i <= length(read); ++i) { // TODO cache length()
            cum_matches_forward[i] = read[i - 1] == gene[i - 1];
            cum_matches_backward[length(read) - i] = read[seqan::length(read) - i] == gene[seqan::length(gene) - i];
        }

        // Compute cumsums TODO Join loops
        for (size_t i = 1; i <= length(read); ++i) { // TODO cache length()
            cum_matches_forward[i] += cum_matches_forward[i - 1];
            cum_matches_backward[seqan::length(read) - i] += cum_matches_backward[seqan::length(read) - i + 1];
        }

        std::vector<int> sum(seqan::length(read) + 1);
        for (size_t i = 0; i <= seqan::length(read); ++i) {
            sum[i] = cum_matches_forward[i] + cum_matches_backward[i];
        }

        return int(std::max_element(sum.cbegin(), sum.cend()) - sum.cbegin());
    }

    using TAlign = seqan::Align<seqan::Dna5String, seqan::ArrayGaps>;     // align type
    TAlign path2seqanAlignment(const AlignmentPath &path,
                               const seqan::Dna5String &read,
                               const seqan::Dna5String &gene);

    template<class T, typename Tf>
    bool is_topologically_sorted(const T &combined, const Tf &has_edge) {
        for (size_t i = 0; i < combined.size(); ++i) {
            for (size_t j = i; j < combined.size(); ++j) {
                if (has_edge(combined[j], combined[i])) {
                    return false;
                }
            }
        }

        return true;
    }
}