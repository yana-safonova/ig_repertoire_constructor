#include "pairwise_block_aligner.hpp"

namespace algorithms {


    std::vector<PairwiseBlockAlignment> PairwiseBlockAligner::query(const seqan::Dna5String &read, size_t limit,
                                                                    size_t start, size_t finish) const {
        auto result = query_unordered(read, start, finish);
        limit = std::min(limit, result.size());
        if (limit == 0) {
            return {};
        }
        using ctuple_type = decltype(*result.cbegin());

        auto score_function = [](const ctuple_type &a) { return a.int_score; };
        auto comp = [&score_function](const ctuple_type &a,
                                      const ctuple_type &b) -> bool { return score_function(a) > score_function(b); };
        // Return top <limit> positions
        std::nth_element(result.begin(), result.begin() + limit, result.end(), comp);
        result.resize(std::min(result.size(), limit));
        std::sort(result.begin(), result.end(), comp);
        return result;
    }

    PairwiseBlockAlignment PairwiseBlockAligner::make_align(const std::vector<Match> &combined,
                                                            const seqan::Dna5String &read,
                                                            const seqan::Dna5String &query,
                                                            size_t needle_index) const {
            auto has_edge = [this](const Match &a, const Match &b) -> bool {
                int read_gap = b.read_pos - a.read_pos;
                int needle_gap = b.needle_pos - a.needle_pos;
                int gap = read_gap - needle_gap;

                if (gap > scoring.max_local_insertions || -gap > scoring.max_local_deletions) return false;

                // Crossing check
                if (a.needle_pos >= b.needle_pos || a.read_pos >= b.read_pos) return false;

                return true;
            };

            auto vertex_weight = [this](const Match &m) -> double {
                return m.length * scoring.match_reward;
            };

            auto edge_weight = [this](const Match &a, const Match &b) -> double {
                int read_gap = b.read_pos - a.read_pos;
                int needle_gap = b.needle_pos - a.needle_pos;
                int gap = read_gap - needle_gap;

                return -Match::overlap(a, b) -
                       ((gap) ? (scoring.gap_opening_cost + std::abs(gap) * scoring.gap_extention_cost) : 0);
            };

            auto _ = weighted_longest_path_in_DAG(combined, has_edge, edge_weight, vertex_weight);
            return PairwiseBlockAlignment::path2Alignment(_.first, read, query, needle_index, _.second);
    }

    bool PairwiseBlockAligner::check_alignment(const PairwiseBlockAlignment &align) const {
            const auto &path = align.path; // FIXME
            if (std::abs(path.global_gap()) >
                scoring.max_global_gap) { // TODO split into 2 args (ins/dels) positive gap is deletion here
                    // Omit such match
                    return false;
            }
            if (path.kplus_length() < min_k_coverage) { // TODO Rename kplus_length()
                    // Omit such match
                    return false;
            }
            return true;
    }

    std::vector<PairwiseBlockAlignment> PairwiseBlockAligner::query_unordered(const seqan::Dna5String &read,
                                                                              size_t start, size_t finish) const {
            std::vector<std::vector<KmerMatch>> needle2matches(queries.size());

            // if (length(read) < K) {
            //     // Return empty result
            //     return {  };
            // }

            finish = std::min(finish, length(read));
            assert(start <= finish);
            if (finish - start < K) {
                    // Return empty result
                    return {};
            }

            auto hashes = polyhashes(read, K);
            for (size_t j = start; j < finish - K + 1; ++j) { // Scan k-mers on given interval
                    auto kmer = hashes[j];
                    auto it = kmer2needle.find(kmer);

                    if (it == kmer2needle.cend()) {
                            continue;
                    }

                    for (const auto &p : it->second) {
                            size_t needle_index = p.needle_index;
                            size_t kmer_pos_in_read = j;
                            size_t kmer_pos_in_needle = p.position;

                            needle2matches[needle_index].push_back(
                                    {static_cast<int>(kmer_pos_in_needle), static_cast<int>(kmer_pos_in_read)});
                    }
            }

            std::vector<PairwiseBlockAlignment> result;

            for (size_t needle_index = 0; needle_index < queries.size(); ++needle_index) {
                    auto &matches = needle2matches[needle_index];

                    if (matches.empty()) continue;

                    std::vector<Match> combined = combine_sequential_kmer_matches(matches, K);
                    std::sort(combined.begin(), combined.end(),
                              [](const Match &a, const Match &b) -> bool { return a.needle_pos < b.needle_pos; });
                    assert(combined.size() > 0);

                    PairwiseBlockAlignment align = make_align(combined, read, queries[needle_index], needle_index);
                    if (check_alignment(align)) {
                            result.push_back(std::move(align));
                    }
            }

            return result;
    }
}