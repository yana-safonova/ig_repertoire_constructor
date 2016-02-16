#pragma once

#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <memory>
#include <unordered_map>
#include <fstream>
#include <boost/algorithm/string/replace.hpp>
#include <mutex>
#include <chrono>
#include <stdexcept>

using std::vector;
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::map;
using std::make_pair;

#include <boost/program_options.hpp>
#include "fast_ig_tools.hpp"
using path::make_dirs;

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::CharString;
using seqan::length;
using seqan::SeqFileIn;
using seqan::SeqFileOut;

#include "ig_matcher.hpp"

namespace fast_ig_tools {

// input vector(pair(needle_pos, read_pos))
// output vector(tuple(needle_pos, read_pos, length))
struct KmerMatch {
    int needle_pos;
    int read_pos;

    int shift() const {
        return read_pos - needle_pos;
    }

    bool operator<(const KmerMatch &m) const {
        return (shift() < m.shift()) || (shift() == m.shift() && read_pos < m.read_pos);
    }
};


// struct Match : public KmerMatch {
struct Match {
    int needle_pos;
    int read_pos;
    size_t length;

    static int overlap(const Match &a,
                       const Match &b) {
        return std::max<int>(std::max<int>(a.length - (b.needle_pos - a.needle_pos),
                                           a.length - (b.read_pos - a.read_pos)),
                             0);
    }
};


int path_coverage_length(const std::vector<Match> &path) {
    int result = 0;

    for (const auto &match : path) {
        result += match.length;
    }

    return result;
}


std::vector<Match> combine_sequential_kmer_matches(std::vector<KmerMatch> &matches,
                                                   size_t K) {
    std::sort(matches.begin(), matches.end());

    std::vector<Match> res;
    res.reserve(matches.size()); // TODO Is it really necessary?

    if (matches.size() == 0) {
        return res;
    }

    Match cur = { matches[0].needle_pos, matches[0].read_pos, K }; // start first match
    for (size_t i = 1; i < matches.size(); ++i) {
        if (matches[i].needle_pos == matches[i-1].needle_pos + 1 && matches[i].read_pos == matches[i-1].read_pos + 1) { // extend current match
            cur.length += 1;
        } else { // save match and start new one
            res.push_back(cur);
            cur = { matches[i].needle_pos, matches[i].read_pos, K };
        }
    }
    res.push_back(cur); // save last match

    return res;
}


std::string visualize_matches(const std::vector<Match> &matches,
                              int needle_length, int read_length) {
    // Draw fancy alignment
    // (read/needle)
    std::stringstream ss;

    assert(std::is_sorted(matches.cbegin(), matches.cend(),
                          [](const Match &a, const Match &b) -> bool { return a.needle_pos < b.needle_pos; }));
    assert(std::is_sorted(matches.cbegin(), matches.cend(),
                          [](const Match &a, const Match &b) -> bool { return a.read_pos < b.read_pos; }));

    ss << bformat("{%d}") % std::max(matches[0].needle_pos - matches[0].read_pos, 0);
    ss << bformat("(%d)") % std::min(matches[0].needle_pos, matches[0].read_pos);
    for (size_t i = 0; i < matches.size() - 1; ++i) {
        int read_gap = matches[i+1].read_pos - matches[i].read_pos - matches[i].length;
        int needle_gap = matches[i+1].needle_pos - matches[i].needle_pos - matches[i].length;
        unsigned current_match_len = matches[i].length;

        if (needle_gap >=0 || read_gap >= 0) {
            if (needle_gap != read_gap) {
                ss << bformat("%d(%d%+d)") % current_match_len % needle_gap % (read_gap - needle_gap);
            } else {
                ss << bformat("%2%(%1%)") % read_gap % current_match_len;
            }
        }
    }

    const auto &last_match = matches[matches.size() - 1];
    ss << bformat("%1%") % last_match.length;
    ss << bformat("(%d)") % std::min(needle_length - last_match.needle_pos - last_match.length,
                                     read_length - last_match.read_pos - last_match.length);
    ss << bformat("{%d}") % std::max((needle_length - last_match.needle_pos) - (read_length - last_match.read_pos), 0);

    return ss.str();
}


class BlockAligner {
    struct PositionInDB {
        size_t needle_index;
        size_t position;
    };

    struct Alignment {
        int kp_coverage;
        std::vector<Match> path;
        int start, finish;
        size_t needle_index;
        int overlap_length;
        double score;

        bool operator< (const Alignment& b) const {
            return this->kp_coverage < b.kp_coverage;
        }


        size_t first_match_read_pos() const {
            return path[0].read_pos;
        }

        size_t first_match_needle_pos() const {
            return path[0].needle_pos;
        }

        size_t last_match_needle_pos() const {
            const auto &last = path[path.size() - 1];
            return last.needle_pos + last.length;
        }

        size_t last_match_read_pos() const {
            const auto &last = path[path.size() - 1];
            return last.read_pos + last.length;
        }
    };

public:
    BlockAligner(const vector<Dna5String> &queries,
                 size_t K,
                 int max_global_gap, int left_uncoverage_limit, int right_uncoverage_limit,
                 int max_local_insertions,
                 int max_local_deletions,
                 int min_k_coverage) : max_local_insertions{max_local_insertions},
                                       max_local_deletions{max_local_deletions},
                                       min_k_coverage{min_k_coverage},
                                       K{K},
                                       left_uncoverage_limit{left_uncoverage_limit},
                                       right_uncoverage_limit{right_uncoverage_limit},
                                       max_global_gap{max_global_gap} {
            this->queries = queries;

            for (size_t j = 0; j < this->queries.size(); ++j) {
                auto hashes = polyhashes(queries[j], K);
                for (size_t start = 0; start + K <= length(queries[j]); ++start) {
                    kmer2needle[hashes[start]].push_back( { j, start } );
                }
            }
        }

    std::vector<Alignment> query(const Dna5String &read, size_t limit) const {
        auto result = query_unordered(read);

        using ctuple_type = decltype(*result.cbegin());

        auto score_function = [](const ctuple_type &a) { return a.kp_coverage; };
        auto comp = [&score_function](const ctuple_type &a,
                                      const ctuple_type &b) -> bool { return score_function(a) > score_function(b); };

        // Return top <limit> positions
        std::nth_element(result.begin(), result.begin() + limit, result.end(), comp);
        result.resize(std::min(result.size(), limit));
        std::sort(result.begin(), result.end(), comp);

        return result;
    }

private:
    std::vector<std::vector<KmerMatch>> needle2matches(const Dna5String &read) const {
        std::vector<std::vector<KmerMatch>> result(queries.size());

        if (length(read) < K) {
            // Return empty result
            return result;
        }

        auto hashes = polyhashes(read, K);
        for (size_t j = 0; j < hashes.size(); ++j) {
            auto kmer = hashes[j];
            auto it = kmer2needle.find(kmer);

            if (it == kmer2needle.cend()) {
                continue;
            }

            for (const auto &p : it->second) {
                size_t needle_index = p.needle_index;
                size_t kmer_pos_in_read = j;
                size_t kmer_pos_in_needle = p.position;
                int shift = static_cast<int>(kmer_pos_in_read) - static_cast<int>(kmer_pos_in_needle);

                // We make these limits less strict because of possibility of indels
                int shift_min = -left_uncoverage_limit - max_global_gap;
                int shift_max = static_cast<int>(length(read)) - static_cast<int>(length(queries[needle_index])) + right_uncoverage_limit + max_global_gap;

                if (shift >= shift_min && shift <= shift_max) {
                    result[needle_index].push_back( { static_cast<int>(kmer_pos_in_needle), static_cast<int>(kmer_pos_in_read) } );
                }
            }
        }

        return result;
    }

    std::vector<Alignment> query_unordered(const Dna5String &read) const {
        std::vector<Alignment> result;
        auto n2m = this->needle2matches(read);

        for (size_t needle_index = 0; needle_index < queries.size(); ++needle_index) {
            auto &matches = n2m[needle_index];

            if (matches.empty()) continue;

            std::vector<Match> combined = combine_sequential_kmer_matches(matches, K);

            assert(combined.size() > 0);

            std::sort(combined.begin(), combined.end(),
                      [](const Match &a, const Match &b) -> bool { return a.needle_pos < b.needle_pos; });

            std::vector<double> values(combined.size(), 0.);
            std::vector<size_t> next(combined.size());
            std::iota(next.begin(), next.end(), 0);

            auto has_edge = [&](const Match &a, const Match &b) -> bool {
                int read_gap = b.read_pos - a.read_pos;
                int needle_gap = b.needle_pos - a.needle_pos;
                int insert_size = read_gap - needle_gap;

                if (insert_size > max_local_insertions || -insert_size > max_local_deletions) return false;

                // Crossing check
                if (a.needle_pos >= b.needle_pos || a.read_pos >= b.read_pos) return false;

                // Obsolete overlap checking
                // return std::min(b.needle_pos - a.needle_pos, b.read_pos - a.read_pos) >= int(a.length); // signed/unsigned comparisson
                return true;
            };

            for (size_t i = combined.size() - 1; i + 1 > 0; --i) {
                values[i] = combined[i].length;

                for (size_t j = i + 1; j < combined.size(); ++j) {
                    if (has_edge(combined[i], combined[j])) {
                        double new_val = combined[i].length + values[j] - Match::overlap(combined[i], combined[j]);
                        if (new_val > values[i]) {
                            next[i] = j;
                            values[i] = new_val;
                        }
                    }
                }
            }

            std::vector<Match> path;
            path.reserve(combined.size());

            size_t maxi = std::max_element(values.cbegin(), values.cend()) - values.cbegin();

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
                path[i].length -= Match::overlap(path[i], path[i+1]);
            }

            int coverage_length = path_coverage_length(path);

            assert(std::is_sorted(path.cbegin(), path.cend(), has_edge));

            // Just use the most left and most right matches
            int left_shift = path[0].read_pos - path[0].needle_pos;
            int right_shift = path[path.size() - 1].read_pos - path[path.size() - 1].needle_pos;

            if (std::abs(left_shift - right_shift) > max_global_gap) {
                // Omit such match
                continue;
            }

            if (coverage_length < min_k_coverage) {
                // Omit such match
                continue;
            }

            int start = left_shift;
            int finish = right_shift + int(length(queries[needle_index]));

            int shift_min = -left_uncoverage_limit;
            int shift_max = int(length(read)) - int(length(queries[needle_index])) + right_uncoverage_limit;

            if (left_shift < shift_min || right_shift > shift_max) {
                // Omit candidates with unproper final shift
                // Maybe we should make these limits less strict because of possibility of indels on edges?
                continue;
            }

            int over_start = std::max(0, start);
            int over_finish = std::min(right_shift + length(queries[needle_index]), length(read));
            int read_overlap_length = over_finish - over_start; // read overlap
            int needle_overlap_length = read_overlap_length + left_shift - right_shift;

            Alignment align;
            align.kp_coverage = coverage_length;
            align.score = static_cast<double>(coverage_length) / static_cast<double>(length(queries[needle_index]));
            align.path = std::move(path);
            align.start = start;
            align.finish = finish;
            align.needle_index = needle_index;
            align.overlap_length = needle_overlap_length;

            result.push_back(std::move(align));
        }

        return result;
    }

    vector<Dna5String> queries;
    int max_local_insertions;
    int max_local_deletions;
    int min_k_coverage;
    size_t K;
    int left_uncoverage_limit, right_uncoverage_limit;
    int max_global_gap;
    std::unordered_map<size_t, std::vector<PositionInDB>> kmer2needle;
};

} // namespace fast_ig_tools

// vim: ts=4:sw=4
