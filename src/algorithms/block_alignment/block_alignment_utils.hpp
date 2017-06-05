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

    template<typename Tf1, typename Tf2, typename Tf3>
    std::pair<AlignmentPath, int> weighted_longest_path_in_DAG(const std::vector<Match> &combined,
                                                               const Tf1 &has_edge,
                                                               const Tf2 &edge_weight,
                                                               const Tf3 &vertex_weight) {
        VERIFY(combined.size() > 0);
        VERIFY(std::is_sorted(combined.cbegin(), combined.cend(), Match::less_subject_pos));
        // Vertices should be topologically sorted
        VERIFY(is_topologically_sorted(combined, has_edge));

        std::vector<double> values(combined.size(), 0.);
        std::vector<size_t> next(combined.size());
        std::iota(next.begin(), next.end(), 0);

        for (size_t i = combined.size() - 1; i + 1 > 0; --i) {
            values[i] = vertex_weight(combined[i]);

            for (size_t j = i + 1; j < combined.size(); ++j) {
                // Check topologically order
                // TODO remove one of these toposort checkings
                assert(!has_edge(combined[j], combined[i]));

                if (has_edge(combined[i], combined[j])) {
                    double new_val = vertex_weight(combined[i]) + values[j] + edge_weight(combined[i], combined[j]);
                    if (new_val > values[i]) {
                        next[i] = j;
                        values[i] = new_val;
                    }
                }
            }
        }

        AlignmentPath path;
        path.reserve(combined.size());

        size_t maxi = size_t(std::max_element(values.cbegin(), values.cend()) - values.cbegin());
        // Sasha, is it ok that score is integer here? Looks like a potential error
        int score = int(values[maxi]);

        while (true) {
            path.push_back(combined[maxi]);
            if (next[maxi] == maxi) {
                break;
            } else {
                maxi = next[maxi];
            }
        }

        // Fix overlaps (truncate tail of left match)
        for (size_t i = 0; i < path.size() - 1; ++i) {
            path[i].length -= Match::overlap(path[i], path[i + 1]);
        }

        // Path should be correct, all edges should be
        VERIFY(std::is_sorted(path.cbegin(), path.cend(), has_edge));

        return { path, score };
    }
}