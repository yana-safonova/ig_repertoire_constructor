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
    protected:
    	const SubjectQueryKmerIndex<SubjectDatabase, StringType> &kmer_index_;
        KmerIndexHelper<SubjectDatabase, StringType> &kmer_index_helper_;
        const BlockAlignmentScoringScheme scoring_;
        const BlockAlignerParams params_;

    public:
        // The following 3 methods were lambda fields. Had to refactor due to gcc 4.8 bugs.
        bool has_edge(const Match &a, const Match &b) const {
            int read_gap = b.read_pos - a.read_pos;
            int needle_gap = b.subject_pos - a.subject_pos;
            int gap = read_gap - needle_gap;
            if (gap > scoring_.max_local_insertions || -gap > scoring_.max_local_deletions) return false;
            // Crossing check
            if (a.subject_pos >= b.subject_pos || a.read_pos >= b.read_pos) return false;
            return true;
        }

        int edge_weight(const Match &a, const Match &b) const {
            int read_gap = b.read_pos - a.read_pos;
            int needle_gap = b.subject_pos - a.subject_pos;
            int gap = read_gap - needle_gap;
            int mmatch = std::min(b.read_pos - a.read_pos - int(a.length),
                                  b.subject_pos - a.subject_pos - int(a.length));
            mmatch = std::max(0, mmatch);
            return - Match::overlap(a, b)
                   - ((gap) ? (scoring_.gap_opening_cost + std::abs(gap) * scoring_.gap_extention_cost) : 0)
                   - ((mmatch) ? scoring_.mismatch_opening_cost + mmatch * scoring_.mismatch_extention_cost : 0);
        }

        double vertex_weight(const Match &m) const {
            return double(m.length) * double(scoring_.match_reward);
        }

        PairwiseBlockAligner(const SubjectQueryKmerIndex<SubjectDatabase, StringType> &kmer_index,
                             KmerIndexHelper<SubjectDatabase, StringType> &kmer_index_helper,
                             BlockAlignmentScoringScheme scoring, BlockAlignerParams params) :
                kmer_index_(kmer_index),
                kmer_index_helper_(kmer_index_helper),
                scoring_(scoring),
                params_(params) { }

        bool CheckAlignment(const PairwiseBlockAlignment &alignment) const {
            // TODO split into 2 args (ins/dels) positive gap is deletion here
            if (std::abs(alignment.path.global_gap()) > scoring_.max_global_gap)
                return false; // Omit such match
            return alignment.path.kplus_length() >= params_.min_kmer_coverage;
        }

        PairwiseBlockAlignment pairwiseBlockAlignment(AlignmentPath &path,
        									size_t subject_index,
        									const StringType &query,
                                            int score) const {
            
            return PairwiseBlockAlignment(path,
                                          kmer_index_helper_.GetStringLength(
                                                  kmer_index_helper_.GetDbRecordByIndex(subject_index)),
                                          kmer_index_helper_.GetStringLength(query),
                                          score);
        }

        virtual BlockAlignmentHits<SubjectDatabase> Align(const StringType &query) = 0;

        virtual ~PairwiseBlockAligner() = default;
    };
    
    template<typename SubjectDatabase, typename StringType>
    class QuadraticDAGPairwiseBlockAligner : public PairwiseBlockAligner<SubjectDatabase, StringType> {
    public:
        QuadraticDAGPairwiseBlockAligner(const SubjectQueryKmerIndex<SubjectDatabase, StringType> &kmer_index,
                             KmerIndexHelper<SubjectDatabase, StringType> &kmer_index_helper,
                             BlockAlignmentScoringScheme scoring, 
                             BlockAlignerParams params) :
        		PairwiseBlockAligner<SubjectDatabase, StringType>(kmer_index, kmer_index_helper, scoring, params) {}

        BlockAlignmentHits<SubjectDatabase> Align(const StringType &query) {
            auto result = QueryUnordered(query);
            result.SelectTopRecords(this->params_.max_candidates);
            return result;
        }

    private:

        using PairwiseBlockAligner<SubjectDatabase, StringType>::has_edge;
        using PairwiseBlockAligner<SubjectDatabase, StringType>::edge_weight;
        using PairwiseBlockAligner<SubjectDatabase, StringType>::vertex_weight;

        PairwiseBlockAlignment MakeAlignment(const std::vector<Match> &combined,
                                          const StringType &query,
                                          size_t subject_index) const {
            
            auto longest_path = weighted_longest_path_in_DAG(combined);
            return this->pairwiseBlockAlignment(longest_path.first,
            							subject_index,
            							query,
                                        longest_path.second);
        }

        BlockAlignmentHits<SubjectDatabase> QueryUnordered(const StringType &query) {
            SubjectKmerMatches subj_matches = this->kmer_index_.GetSubjectKmerMatchesForQuery(query);
            BlockAlignmentHits<SubjectDatabase> result(this->kmer_index_.Db());
            for(size_t i = 0; i < subj_matches.size(); i++) {
                auto &matches = subj_matches[i];
                if(matches.empty())
                    continue;
                std::vector<Match> combined = combine_sequential_kmer_matches(matches, this->kmer_index_.k());
                std::sort(combined.begin(), combined.end(),
                          [](const Match &a, const Match &b) -> bool { return a.subject_pos < b.subject_pos; });
                VERIFY(combined.size() > 0);
                PairwiseBlockAlignment align = MakeAlignment(combined, query, i);
                if (this->CheckAlignment(align)) {
                    result.Add(std::move(align), i);
                }
            }
            return result;
        }

        std::pair<AlignmentPath, int> weighted_longest_path_in_DAG(const std::vector<Match> &combined) const {
            VERIFY(!combined.empty());
            VERIFY(std::is_sorted(combined.cbegin(), combined.cend(), Match::less_subject_pos));
            // Vertices should be topologically sorted
            VERIFY(is_topologically_sorted(combined, [this](const Match &a, const Match &b) -> bool { return has_edge(a, b); }));

            std::vector<double> values(combined.size(), 0.);
            std::vector<size_t> next(combined.size());
            std::iota(next.begin(), next.end(), 0);

            for (size_t i = combined.size() - 1; i + 1 > 0; --i) {
                values[i] = vertex_weight(combined[i]);

                for (size_t j = i + 1; j < combined.size(); ++j) {
                    // Check topologically order
                    // TODO remove one of these toposort checkings
                    VERIFY(!has_edge(combined[j], combined[i]));

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
            VERIFY(std::is_sorted(path.cbegin(), path.cend(), [this](const Match &a, const Match &b) -> bool { return has_edge(a, b); }));

            return { path, score };
        }
    };

    template<typename SubjectDatabase, typename StringType>
    class LisPairwiseBlockAligner : public PairwiseBlockAligner<SubjectDatabase, StringType> {

    public:
        LisPairwiseBlockAligner(const SubjectQueryKmerIndex<SubjectDatabase, StringType> &kmer_index,
                             KmerIndexHelper<SubjectDatabase, StringType> &kmer_index_helper,
                             BlockAlignmentScoringScheme scoring, BlockAlignerParams params) :
                PairwiseBlockAligner<SubjectDatabase, StringType>(kmer_index, kmer_index_helper, scoring, params) {}

        BlockAlignmentHits<SubjectDatabase> Align(const StringType &query) {
            auto result = QueryUnordered(query);
            result.SelectTopRecords(this->params_.max_candidates);
            return result;
        }

    private:
    	BlockAlignmentHits<SubjectDatabase> QueryUnordered(const StringType &query) {
            SubjectKmerMatches subj_matches = this->kmer_index_.GetSubjectKmerMatchesForQuery(query);
            BlockAlignmentHits<SubjectDatabase> result(this->kmer_index_.Db());
            for(size_t i = 0; i < subj_matches.size(); i++) {
                auto &matches = subj_matches[i];
                if(matches.empty())
                    continue;
                std::sort(matches.begin(), matches.end(),
                          [](const KmerMatch &a, const KmerMatch &b) -> bool { return a.needle_pos < b.needle_pos; });
                matches.resize(static_cast<size_t>(
                        std::unique(matches.begin(), matches.end(),
                                    [](const KmerMatch &a, const KmerMatch &b) -> bool { return a.needle_pos == b.needle_pos; }
                                    ) - matches.begin()
                               ));

                PairwiseBlockAlignment align = MakeAlignment(matches, query, i);
                if (this->CheckAlignment(align)) {
                    result.Add(std::move(align), i);
                }
            }
            return result;
        }

        std::vector<KmerMatch> findLIS(const std::vector<KmerMatch> &matches) const {

            std::vector<int> pos_before(matches.size() + 1, -1);
            std::vector<int> pos(matches.size() + 1, -1);
            std::vector<int> min_value(matches.size() + 1, std::numeric_limits<int>::max());
            min_value[0] = std::numeric_limits<int>::min();
            int max_len = 0;
            for (size_t i = 0; i < matches.size(); i++) {
                int new_value = matches[i].read_pos;
                int cur_len = static_cast<int>(std::lower_bound(begin(min_value), end(min_value), new_value) - begin(min_value));
                min_value[cur_len] = new_value;
                pos[cur_len] = (int) i;
                pos_before[i] = pos[cur_len - 1];
                max_len = std::max(cur_len, max_len);
            }
            std::vector<KmerMatch> res;
            {
                int i = pos[max_len];
                while (i >= 0) {
                    res.push_back(matches[i]);
                    i = pos_before[i];
                }
            }

            std::reverse(res.begin(), res.end());
            for (size_t i = 1; i < res.size(); i++) {
                VERIFY(res[i - 1].needle_pos <= res[i].needle_pos);
                VERIFY(res[i - 1].read_pos < res[i].read_pos);
            }
            return res;
        }

        PairwiseBlockAlignment MakeAlignment(const std::vector<KmerMatch> &matches,
                                          const StringType &query,
                                          size_t subject_index) const {

            std::vector<KmerMatch> lis = findLIS(matches);
            std::vector<Match> combined = combine_sequential_kmer_matches(lis, this->kmer_index_.k());
            std::sort(begin(combined), end(combined), [](const Match a, const Match b) -> bool {
                return a.read_pos < b.read_pos;
            });
            for (size_t i = 1; i < combined.size(); i++) {
                VERIFY(combined[i - 1].subject_pos < combined[i].subject_pos);
                VERIFY(combined[i - 1].read_pos < combined[i].read_pos);
            }

            AlignmentPath path;
            for (size_t i = 0; i < combined.size(); i++) { 
                path.push_back(combined[i]);
            }
            for (size_t i = 0; i + 1 < path.size(); ++i) {
                path[i].length -= Match::overlap(path[i], path[i + 1]);
            }

            double score = 0;
            for (size_t i = 0; i < combined.size(); i++) { 
                score += this->vertex_weight(combined[i]);
            }
            for (size_t i = 1; i < combined.size(); i++) {
                score += this->edge_weight(combined[i - 1], combined[i]);
            }

            return this->pairwiseBlockAlignment(path,
                                          subject_index,
                                          query,
                                          (int) score);
        }
    };

	template<typename SubjectDatabase, typename StringType>
    class QuadraticDpPairwiseBlockAligner : public PairwiseBlockAligner<SubjectDatabase, StringType> {

   	public:
        QuadraticDpPairwiseBlockAligner(const SubjectQueryKmerIndex<SubjectDatabase, StringType> &kmer_index,
                             KmerIndexHelper<SubjectDatabase, StringType> &kmer_index_helper,
                             BlockAlignmentScoringScheme scoring, BlockAlignerParams params) :
                PairwiseBlockAligner<SubjectDatabase, StringType>(kmer_index, kmer_index_helper, scoring, params) {}

        BlockAlignmentHits<SubjectDatabase> Align(const StringType &query) {
            auto result = QueryUnordered(query);
            result.SelectTopRecords(this->params_.max_candidates);
            return result;
        }

    private:
    	BlockAlignmentHits<SubjectDatabase> QueryUnordered(const StringType &query) {
            SubjectKmerMatches subj_matches = this->kmer_index_.GetSubjectKmerMatchesForQuery(query);
            BlockAlignmentHits<SubjectDatabase> result(this->kmer_index_.Db());
            for(size_t i = 0; i < subj_matches.size(); i++) {
                auto &matches = subj_matches[i];
                if(matches.empty())
                    continue;
                std::vector<Match> combined = combine_sequential_kmer_matches(matches, this->kmer_index_.k());
                assert(combined.size() > 0);
                PairwiseBlockAlignment align = MakeAlignment(combined, query, i);
                if (this->CheckAlignment(align)) {
                    result.Add(std::move(align), i);
                }
            }
            return result;
        }

        PairwiseBlockAlignment MakeAlignment(const std::vector<Match> &combined,
                                          const StringType &query,
                                          size_t subject_index) const {

        	std::vector<int> dp_order(combined.size());
        	for (size_t i = 0; i < dp_order.size(); i++) {
        		dp_order[i] = (int)i;
        	}

        	std::sort(dp_order.begin(), dp_order.end(), [&combined](const int i1, const int i2) -> bool{
        		return combined[i1].subject_pos < combined[i2].subject_pos;
        	});

            std::vector<int> dp(combined.size());
        	std::vector<int> pos_before(combined.size());
        	int mx_dp_pos = 0;

        	for (size_t i = 0; i < combined.size(); i++) {
        		int pos = dp_order[i];
        		int mx = 0;
        		int mx_pos = -1;
        		const Match &b = combined[pos];

        		for (size_t pos2 = static_cast<size_t>(pos + 1); pos2 < combined.size(); pos2++) {
        			const Match &a = combined[pos2];
        			int read_gap = b.read_pos - a.read_pos;
		        	int needle_gap = b.subject_pos - a.subject_pos;
		        	int gap = read_gap - needle_gap;
		        	if (gap > this->scoring_.max_local_insertions || -gap > this->scoring_.max_local_deletions) {
		        		break;
		        	}
		

                	if (this->has_edge(a, b)) {
                		int cost = this->edge_weight(a, b);
                		if (cost > mx) {
                			mx = cost;
                			mx_pos = (int) pos2;
                		}
        			}
        		}

        		for (size_t pos2 = static_cast<size_t>(pos - 1); pos2 + 1 > 0; pos2--) {
        			const Match &a = combined[pos2];
        			int read_gap = b.read_pos - a.read_pos;
		        	int needle_gap = b.subject_pos - a.subject_pos;
		        	int gap = read_gap - needle_gap;
		        	if (gap > this->scoring_.max_local_insertions || -gap > this->scoring_.max_local_deletions) {
		        		break;
		        	}

                	if (this->has_edge(a, b)) {
                		int cost = this->edge_weight(a, b) + dp[pos2];
                		if (cost > mx) {
                			mx = cost;
                			mx_pos = (int) pos2;
                		}
        			}
        		}

        		dp[pos] = (int)this->vertex_weight(b) + mx;
        		pos_before[pos] = mx_pos;

        		if (dp[pos] > dp[mx_dp_pos]) {
        			mx_dp_pos = pos;
        		}
        	}

        	int pos = mx_dp_pos;
        	int score = dp[mx_dp_pos];
        	AlignmentPath path;
        	while (pos != -1) {
        		path.push_back(combined[pos]);
        		pos = pos_before[pos];
        	}
            std::reverse(path.begin(), path.end());
            for (size_t i = 0; i + 1 < path.size(); ++i) {
                path[i].length -= Match::overlap(path[i], path[i + 1]);
            }

            VERIFY(std::is_sorted(path.cbegin(), path.cend(), [this](const Match &a, const Match &b) -> bool { return this->has_edge(a, b); }));

            return this->pairwiseBlockAlignment(path,
                                          subject_index,
                                          query,
                                          score);
        }
    };
}
