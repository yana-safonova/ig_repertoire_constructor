#pragma once

#include "pairwise_block_alignment.hpp"
#include "../hashes/subject_query_kmer_index.hpp"

namespace algorithms {
    struct BlockAlignmentScoringScheme {
        int max_global_gap;
        int max_local_deletions;
        int max_local_insertions;
        int gap_opening_cost;
        int gap_extention_cost;
        int match_reward;
        int mismatch_extention_cost;
        int mismatch_opening_cost;
    };

    struct BlockAlignerParams {
        size_t min_kmer_coverage;
        size_t max_candidates;

        BlockAlignerParams(size_t min_kmer_coverage, size_t max_candidates) :
                min_kmer_coverage(min_kmer_coverage),
                max_candidates(max_candidates) { }
    };

    std::vector<Match> combine_sequential_kmer_matches(std::vector<KmerMatch> &matches,
                                                       size_t K);

    template<typename SubjectDatabase, typename StringType>
    class PairwiseBlockAligner {
        const SubjectQueryKmerIndex<SubjectDatabase, StringType> &kmer_index_;
        KmerIndexHelper<SubjectDatabase, StringType> &kmer_index_helper_;
        const BlockAlignmentScoringScheme scoring_;
        const BlockAlignerParams params_;

        PairwiseBlockAlignment MakeAlignment(const std::vector<Match> &combined,
                                          const StringType &query,
                                          size_t subject_index) const {
            auto has_edge = [this](const Match &a, const Match &b) -> bool {
                int read_gap = b.read_pos - a.read_pos;
                int needle_gap = b.subject_pos - a.subject_pos;
                int gap = read_gap - needle_gap;
                if (gap > scoring_.max_local_insertions || -gap > scoring_.max_local_deletions) return false;
                // Crossing check
                if (a.subject_pos >= b.subject_pos || a.read_pos >= b.read_pos) return false;
                return true;
            };

            auto vertex_weight = [this](const Match &m) -> double {
                return double(m.length) * double(scoring_.match_reward);
            };

            auto edge_weight = [this](const Match &a, const Match &b) -> double {
                int read_gap = b.read_pos - a.read_pos;
                int needle_gap = b.subject_pos - a.subject_pos;
                int gap = read_gap - needle_gap;
                int mmatch = std::min(b.read_pos - a.read_pos - int(a.length),
                                      b.subject_pos - a.subject_pos - int(a.length));
                mmatch = std::max(0, mmatch);
                return - Match::overlap(a, b)
                       - ((gap) ? (scoring_.gap_opening_cost + std::abs(gap) * scoring_.gap_extention_cost) : 0)
                       - ((mmatch) ? scoring_.mismatch_opening_cost + mmatch * scoring_.mismatch_extention_cost : 0);
            };

            auto longest_path = weighted_longest_path_in_DAG(combined, has_edge, edge_weight, vertex_weight);
            return PairwiseBlockAlignment(longest_path.first,
                                          kmer_index_helper_.GetStringLength(
                                                  kmer_index_helper_.GetDbRecordByIndex(subject_index)),
                                          kmer_index_helper_.GetStringLength(query),
                                          longest_path.second);
        }

        bool CheckAlignment(const PairwiseBlockAlignment &alignment) const {
            // TODO split into 2 args (ins/dels) positive gap is deletion here
            if (std::abs(alignment.path.global_gap()) > scoring_.max_global_gap)
                return false; // Omit such match
            return alignment.path.kplus_length() >= params_.min_kmer_coverage;
        }

        BlockAlignmentHits<SubjectDatabase> QueryUnordered(const StringType &query) {
            SubjectKmerMatches subj_matches = kmer_index_.GetSubjectKmerMatchesForQuery(query);
            BlockAlignmentHits<SubjectDatabase> result(kmer_index_.Db());
            for(size_t i = 0; i < subj_matches.size(); i++) {
                auto &matches = subj_matches[i];
                if(matches.empty())
                    continue;
                std::vector<Match> combined = combine_sequential_kmer_matches(matches, kmer_index_.k());
                std::sort(combined.begin(), combined.end(),
                          [](const Match &a, const Match &b) -> bool { return a.subject_pos < b.subject_pos; });
                assert(combined.size() > 0);
                PairwiseBlockAlignment align = MakeAlignment(combined, query, i);
                if (CheckAlignment(align)) {
                    result.Add(std::move(align), i);
                }
            }
            return result;
        }

    public:
        PairwiseBlockAligner(const SubjectQueryKmerIndex<SubjectDatabase, StringType> &kmer_index,
                             KmerIndexHelper<SubjectDatabase, StringType> &kmer_index_helper,
                             BlockAlignmentScoringScheme scoring, BlockAlignerParams params) :
                kmer_index_(kmer_index),
                kmer_index_helper_(kmer_index_helper),
                scoring_(scoring),
                params_(params) { }

        BlockAlignmentHits<SubjectDatabase> Align(const StringType &query) {
            auto result = QueryUnordered(query);
            result.SelectTopRecords(params_.max_candidates);
            return result;
        }
    };
}
