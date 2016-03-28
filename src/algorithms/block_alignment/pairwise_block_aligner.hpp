#pragma once

#include "pairwise_block_alignment.hpp"
#include "../hashes/polyhashes.hpp"

namespace algorithms {
    class PairwiseBlockAligner {
        struct PositionInDB {
            size_t needle_index;
            size_t position;
        };

    public:
        struct ScoringScheme {
            int max_global_gap = 24;
            int max_local_deletions = 12;
            int max_local_insertions = 12;

            int gap_opening_cost = 4;
            int gap_extention_cost = 1;
            int match_reward = 1;
        };


        PairwiseBlockAligner(const std::vector<seqan::Dna5String> &queries,
                     size_t K,
                     const ScoringScheme &scoring,
                     int min_k_coverage) : K{K},
                                           scoring{scoring},
                                           min_k_coverage{min_k_coverage} {
            this->queries = queries;

            for (size_t j = 0; j < this->queries.size(); ++j) {
                auto hashes = polyhashes(queries[j], K);
                for (size_t start = 0; start + K <= length(queries[j]); ++start) {
                    kmer2needle[hashes[start]].push_back({j, start});
                }
            }
        }

        std::vector<PairwiseBlockAlignment> query(const seqan::Dna5String &read,
                                     size_t limit,
                                     size_t start = 0, size_t finish = 10005000) const;

    private:
        PairwiseBlockAlignment make_align(const std::vector<Match> &combined,
                             const seqan::Dna5String &read,
                             const seqan::Dna5String &query,
                             size_t needle_index) const;

        bool check_alignment(const PairwiseBlockAlignment &align) const;

        std::vector<PairwiseBlockAlignment> query_unordered(const seqan::Dna5String &read,
                                               size_t start = 0, size_t finish = 10005000) const;

        std::vector<seqan::Dna5String> queries;
        size_t K;
        ScoringScheme scoring;
        int min_k_coverage;
        std::unordered_map<size_t, std::vector<PositionInDB>> kmer2needle;
    };
}