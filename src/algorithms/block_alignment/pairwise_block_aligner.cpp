#include "pairwise_block_aligner.hpp"

namespace algorithms {
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
            if (matches[i].needle_pos == matches[i-1].needle_pos + 1 &&
                matches[i].read_pos == matches[i-1].read_pos + 1) { // extend current match
                cur.length += 1;
            } else { // save match and start new one
                res.push_back(cur);
                cur = { matches[i].needle_pos, matches[i].read_pos, K };
            }
        }
        res.push_back(cur); // save last match
        return res;
    }
}