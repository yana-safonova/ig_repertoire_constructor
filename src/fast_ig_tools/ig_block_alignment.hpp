#pragma once

#include <cassert>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <stdexcept>

#include "fast_ig_tools.hpp"

#include <seqan/seq_io.h>
#include <seqan/align.h>

#include "ig_matcher.hpp"

namespace fast_ig_tools {
using seqan::Dna5String;
using seqan::length;

// input vector(pair(needle_pos, read_pos))
// output vector(tuple(needle_pos, read_pos, length))
struct KmerMatch {
    int needle_pos;
    int read_pos;

    int shift() const {
        return read_pos - needle_pos;
    }

    static bool less_shift(const KmerMatch &m1, const KmerMatch &m2) {
        return (m1.shift() < m2.shift()) || (m1.shift() == m2.shift() && m1.read_pos < m2.read_pos);
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

    static bool less_needle_pos(const Match &a, const Match &b) {
        return a.needle_pos < b.needle_pos;
    }

    static bool less_read_pos(const Match &a, const Match &b) {
        return a.read_pos < b.read_pos;
    }
};

class AlignmentPath : public std::vector<Match> {
    using std::vector<Match>::vector;
public:
    int kplus_length() const {
        int result = 0;

        for (const auto &match : *this) {
            result += match.length;
        }

        return result;
    }

    const Match& first() const {
        return (*this)[0];
    }

    const Match& last() const {
        return (*this)[this->size() - 1];
    }

    int left_shift() const {
        return first().read_pos - first().needle_pos;
    }

    int right_shift() const {
        return last().read_pos - last().needle_pos;
    }

    int global_gap() const {
        return left_shift() - right_shift();
    }

    int read_segment_size() const {
        return last().read_pos - first().read_pos + last().length;
    }

    int needle_segment_size() const {
        return last().needle_pos - first().needle_pos + last().length;
    }

    std::string visualize_matches(int needle_length, int read_length) const {
        const std::vector<Match> &matches = *this; // TODO Fixit

        // Draw fancy alignment
        // (read/needle)
        std::stringstream ss;

        assert(std::is_sorted(matches.cbegin(), matches.cend(), Match::less_needle_pos));
        assert(std::is_sorted(matches.cbegin(), matches.cend(), Match::less_read_pos));

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

    bool check_overlaps() const {
        if (size() > 1) {
            for (size_t i = 0; i < size() -1; ++i) {
                if (Match::overlap((*this)[i], (*this)[i + 1])) {
                    return false;
                }
            }
        }

        return true;
    }
};


template<typename Tsequence1, typename Tsequence2>
int find_simple_gap(const Tsequence1 &read, const Tsequence2 &gene) {
    // Positive if gaps in read (|read| < |gene|)
    // ggggggggg
    // rrrrrr
    //    rrrrrr
    assert(length(read) < length(gene));

    std::vector<int> cum_matches_forward(length(read) + 1), cum_matches_backward(length(read) + 1);

    for (size_t i = 1; i <= length(read); ++i) { // TODO cache length()
        cum_matches_forward[i] = read[i - 1] == gene[i - 1];
        cum_matches_backward[length(read) - i] = read[length(read) - i] == gene[length(gene) - i];
    }

    // Compute cumsums TODO Join loops
    for (size_t i = 1; i <= length(read); ++i) { // TODO cache length()
        cum_matches_forward[i] += cum_matches_forward[i - 1];
        cum_matches_backward[length(read) - i] += cum_matches_backward[length(read) - i + 1];
    }

    // for (const auto &_ : cum_matches_forward) {
    //     std::cout << _ << " ";
    // }
    // std::cout << " / ";
    // for (const auto &_ : cum_matches_backward) {
    //     std::cout << _ << " ";
    // }
    // std::cout << std::endl;


    std::vector<int> sum(length(read) + 1);
    for (size_t i = 0; i <= length(read); ++i) {
        sum[i] = cum_matches_forward[i] + cum_matches_backward[i];
    }

    return std::max_element(sum.cbegin(), sum.cend()) - sum.cbegin();
}


using TAlign = seqan::Align<Dna5String, seqan::ArrayGaps>;     // align type
TAlign path2seqanAlignment(const AlignmentPath &path, const Dna5String &read, const Dna5String &gene,
                           size_t clipped_head = 0) {
    assert(path.check_overlaps());
    assert(!path.empty());

    using namespace seqan;
    using TRow =  seqan::Row<TAlign>::Type; // gapped sequence type

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), gene);
    assignSource(row(align, 1), read);

    TRow & row_gene = row(align, 0);
    TRow & row_read = row(align, 1);

    // Clip head
    setBeginPosition(row_read, clipped_head);
    size_t read_len = length(read) - clipped_head;
    size_t gene_len = length(gene);

    // Add finishing gaps (reverse order since insert_gaps works with VIEW position!)
    // TODO Use SeqAn clipping if possible
    int finishing_gap = (gene_len - path.last().needle_pos) - (read_len - path.last().read_pos);
    if (finishing_gap > 0) {
        insertGaps(row_read, read_len, finishing_gap);
    } else if (finishing_gap < 0) {
        // insertGaps(row_gene, gene_len, -finishing_gap);
        // Use clipping:
        setEndPosition(row_read, clipped_head + read_len + finishing_gap);
    } else {
        // Do nothing
    }

    // Add edge gaps if needed. Reverse order!
    if (path.size() > 1) {
        for (size_t i = path.size() - 2; i + 1 > 0; --i) {
            const auto &read_edge = infix(read, path[i].read_pos + path[i].length, path[i + 1].read_pos);
            const auto &gene_edge = infix(gene, path[i].needle_pos + path[i].length, path[i + 1].needle_pos);

            // INFO("EDGE: " << read_edge << " - " << gene_edge);

            if (length(read_edge) < length(gene_edge)) {
                insertGaps(row_read,
                           path[i].read_pos + path[i].length + find_simple_gap(read_edge, gene_edge),
                           length(gene_edge) - length(read_edge));
            } else if (length(gene_edge) < length(read_edge)) {
                insertGaps(row_gene,
                           path[i].needle_pos + path[i].length + find_simple_gap(gene_edge, read_edge),
                           length(read_edge) - length(gene_edge));
            } else {
                // Don't add gaps
            }
        }
    }

    // Add starting gaps (reverse order since insert_gaps works with VIEW position!)
    int starting_gap = path.first().needle_pos - path.first().read_pos;
    if (starting_gap > 0) {
        insertGaps(row_read, 0, starting_gap);
    } else if (starting_gap < 0) {
        // insertGaps(row_gene, 0, -starting_gap);
        // Use clipping:
        setBeginPosition(row_read, clipped_head - starting_gap);
    } else {
        // Do nothing
    }

    return align;
}

// void replace_prefix_for_full_read(TAlign &align, const Dna5String &read) {
//     using namespace seqan;
//     using TRow =  seqan::Row<TAlign>::Type; // gapped sequence type
//
//     TRow & row_gene = row(align, 0);
//     TRow & row_read = row(align, 1);
//
//     size_t old_len = length(source(row_read));
//     assert(old_len <= length(read));
//
//     assignSource(row_read, read);
//     setBeginPosition(row_read, beginPosition(row_read) + length(read) - old_len);
// }


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
AlignmentPath weighted_longest_path_in_DAG(const std::vector<Match> &combined,
                                           const Tf1 &has_edge,
                                           const Tf2 &edge_weight,
                                           const Tf3 &vertex_weight) {
    assert(combined.size() > 0);

    assert(std::is_sorted(combined.cbegin(), combined.cend(), Match::less_needle_pos));

    // Vertices should be topologically sorted
    assert(is_topologically_sorted(combined, has_edge));

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
        path[i].length -= Match::overlap(path[i], path[i + 1]);
    }

    // Path should be correct, all edges should be
    assert(std::is_sorted(path.cbegin(), path.cend(), has_edge));

    return path;
}


std::vector<Match> combine_sequential_kmer_matches(std::vector<KmerMatch> &matches,
                                                   size_t K) {
    std::sort(matches.begin(), matches.end(), KmerMatch::less_shift);

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


class BlockAligner {
    struct PositionInDB {
        size_t needle_index;
        size_t position;
    };

public:
    struct Alignment {
        int kp_coverage;
        AlignmentPath path;
        int start, finish;
        size_t needle_index;
        size_t needle_length;
        int overlap_length;
        double score;

        bool operator< (const Alignment& b) const {
            return this->kp_coverage < b.kp_coverage;
        }

        size_t first_match_read_pos() const {
            return path.first().read_pos;
        }

        size_t first_match_needle_pos() const {
            return path.first().needle_pos;
        }

        size_t last_match_needle_pos() const {
            return path.last().needle_pos + path.last().length;
        }

        size_t last_match_read_pos() const {
            return path.last().read_pos + path.last().length;
        }

        size_t left_half_segment_length() const {
            return last_match_needle_pos();
        }

        size_t right_half_segment_length() const {
            return needle_length - first_match_needle_pos();
        }

        size_t segment_length() const {
            return last_match_needle_pos() - first_match_needle_pos();
        }


        TAlign seqan_alignment(const Dna5String &read,
                               const Dna5String &gene,
                               size_t clipped_head = 0) const {
            return path2seqanAlignment(this->path, read, gene, clipped_head);
        }

        static Alignment path2Alignment(AlignmentPath &path,
                                        const Dna5String &read,
                                        const Dna5String &query,
                                        size_t needle_index) {
            int coverage_length = path.kplus_length();

            int left_shift = path.left_shift();
            int right_shift = path.right_shift();

            int start = left_shift;
            int finish = right_shift + int(length(query));

            int over_start = std::max(0, start);
            int over_finish = std::min(right_shift + length(query), length(read));
            int read_overlap_length = over_finish - over_start; // read overlap
            int needle_overlap_length = read_overlap_length + left_shift - right_shift;

            Alignment align;
            align.kp_coverage = coverage_length;
            align.score = static_cast<double>(coverage_length) / static_cast<double>(length(query));
            align.path = std::move(path);
            align.start = start;
            align.finish = finish;
            align.needle_index = needle_index;
            align.overlap_length = needle_overlap_length;
            align.needle_length = length(query);

            return align;
        }
    };

public:
    BlockAligner(const std::vector<Dna5String> &queries,
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

    std::vector<Alignment> query(const Dna5String &read, size_t limit, size_t start = 0, size_t finish = 10005000) const {
        auto result = query_unordered(read, start, finish);

        // std::cout << "SIZE: " << result.size() << std::endl;

        limit = std::min(limit, result.size());
        if (limit == 0) {
            return {  };
        }

        using ctuple_type = decltype(*result.cbegin());

        auto score_function = [](const ctuple_type &a) { return a.kp_coverage; };
        auto comp = [&score_function](const ctuple_type &a,
                                      const ctuple_type &b) -> bool { return score_function(a) > score_function(b); };

        // Return top <limit> positions
        std::nth_element(result.begin(), result.begin() + limit, result.end(), comp);
        result.resize(std::min(result.size(), limit));
        std::sort(result.begin(), result.end(), comp);

        // const auto &gene = this->queries[result[0].needle_index];
        // std::cout << read << std::endl;
        // std::cout << gene << std::endl;
        // std::cout << result[0].path.visualize_matches(length(gene), length(read)) << std::endl;
        // auto align = path2seqanAlignment(result[0].path, read, gene);
        // // replace_prefix_for_full_read(align, "AAAAAAAAAAAAAAAAAA" + read);
        //
        // Dna5String long_read = "AAAAAAAAA";
        // int clipped_head = length(long_read);
        // long_read += read;
        // auto align2 = path2seqanAlignment(result[0].path, long_read, gene, clipped_head);
        // std::cout << align;
        // std::cout << align2;

        return result;
    }

private:
    Alignment make_align(const std::vector<Match> &combined, const Dna5String &read, const Dna5String &query, size_t needle_index) const {
        // std::sort(combined.begin(), combined.end(),
        //           [](const Match &a, const Match &b) -> bool { return a.needle_pos < b.needle_pos; });
        // assert(combined.size() > 0);
        //
        auto has_edge = [this](const Match &a, const Match &b) -> bool {
            int read_gap = b.read_pos - a.read_pos;
            int needle_gap = b.needle_pos - a.needle_pos;
            int insert_size = read_gap - needle_gap;

            if (insert_size > max_local_insertions || -insert_size > max_local_deletions) return false;

            // Crossing check
            if (a.needle_pos >= b.needle_pos || a.read_pos >= b.read_pos) return false;

            return true;
        };

        auto vertex_weight = [](const Match &m) -> double {
            return m.length;
        };

        auto edge_weight = [&](const Match &a, const Match &b) -> double {
            return -Match::overlap(a, b);
        };

        auto path = weighted_longest_path_in_DAG(combined, has_edge, edge_weight, vertex_weight);

        return Alignment::path2Alignment(path, read, query, needle_index);
    }

    bool check_alignment(const Alignment &align, const Dna5String &read, const Dna5String &query,
                         size_t start, size_t finish) const {
        const auto &path =  align.path; // FIXME

        if (std::abs(path.global_gap()) > max_global_gap) { // TODO split into 2 args (ins/dels) positive gap is deletion here
            // Omit such match
            return false;
        }

        if (path.kplus_length() < min_k_coverage) { // TODO Rename kplus_length()
            // Omit such match
            return false;
        }

        int shift_min = start - left_uncoverage_limit;
        int shift_max = static_cast<int>(finish) - static_cast<int>(length(query)) + right_uncoverage_limit;

        if (path.left_shift() < shift_min || path.right_shift() > shift_max) {
            // Omit candidates with unproper final shift
            // Maybe we should make these limits less strict because of possibility of indels on edges?
            return false;
        }

        return true;
    }

    std::vector<Alignment> query_unordered(const Dna5String &read, size_t start = 0, size_t finish = 10005000) const {
        std::vector<std::vector<KmerMatch>> needle2matches(queries.size());

        // if (length(read) < K) {
        //     // Return empty result
        //     return {  };
        // }

        finish = std::min(finish, length(read));
        assert(start <= finish);
        if (finish - start < K) {
            // Return empty result
            return {  };
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
                int shift = static_cast<int>(kmer_pos_in_read) - static_cast<int>(kmer_pos_in_needle);

                // We make these limits less strict because of possibility of indels
                // int shift_min = -left_uncoverage_limit - max_global_gap;
                // int shift_max = static_cast<int>(length(read)) - static_cast<int>(length(queries[needle_index])) + right_uncoverage_limit + max_global_gap;
                int shift_min = static_cast<int>(start) - left_uncoverage_limit - max_global_gap;
                int shift_max = static_cast<int>(finish) - static_cast<int>(length(queries[needle_index])) + right_uncoverage_limit + max_global_gap;

                if (shift >= shift_min && shift <= shift_max) {
                    needle2matches[needle_index].push_back( { static_cast<int>(kmer_pos_in_needle), static_cast<int>(kmer_pos_in_read) } );
                }
            }
        }

        std::vector<Alignment> result;

        for (size_t needle_index = 0; needle_index < queries.size(); ++needle_index) {
            auto &matches = needle2matches[needle_index];

            if (matches.empty()) continue;

            std::vector<Match> combined = combine_sequential_kmer_matches(matches, K);
            std::sort(combined.begin(), combined.end(),
                      [](const Match &a, const Match &b) -> bool { return a.needle_pos < b.needle_pos; });
            assert(combined.size() > 0);

            Alignment align = make_align(combined, read, queries[needle_index], needle_index);
            if (check_alignment(align, read, queries[needle_index], start, finish)) {
                result.push_back(std::move(align));
            }
        }

        return result;
    }

    std::vector<Dna5String> queries;
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
