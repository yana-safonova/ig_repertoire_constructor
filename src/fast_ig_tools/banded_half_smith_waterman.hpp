// half Smith-Waterman distance

#pragma once

// d("", "") = 0
// d(s, "") = d("", s) = lizard_tail(|s|)
// d(Xs, Yz) = max(match(X, Y) + d(s, z); indel_cost + d(Xs, z); indel_cost + d(s, Yz))
#include <seqan/seq_io.h>
#include <cassert>
#include <vector>
#include <algorithm>


template<typename Ts1, typename Ts2, typename Tf>
int half_hamming(const Ts1 &s1, const Ts2 &s2,
                 int match, int mismatch,
                 const Tf &lizard_tail) {
    auto len1 = length(s1), len2 = length(s2); // Cache lengths because of length(s) computation maybe not O(1) e.g. for char*-strings
    auto len = std::min(len1, len2);

    int res = 0;
    for (size_t i = 0; i < len; ++i) {
        res += (s1[i] == s2[i]) ? match : mismatch;
    }

    res += lizard_tail(std::abs(static_cast<int>(len1) - static_cast<int>(len2)));

    return res;
}


template<typename Ts1, typename Ts2, typename Tf>
int half_sw_banded(const Ts1 &s1, const Ts2 &s2,
                   int match, int mismatch,
                   int indel,
                   const Tf &lizard_tail,
                   int max_indels = 6) {
    using std::min;
    using std::max;
    using seqan::length;
    using std::vector;

    const int INF = 1005000;

    if (max_indels == 0) {
        return half_hamming(s1, s2, match, mismatch, lizard_tail);
    }

    int len1 = static_cast<int>(length(s1));  // Cache lengths because of length(s) computation maybe not O(1) e.g. for char*-strings
    int len2 = static_cast<int>(length(s2));  // Cache lengths because of length(s) computation maybe not O(1) e.g. for char*-strings

    vector<int> base(2*max_indels + 1, -INF);
    for (auto i2 = max(0, len1 - max_indels); i2 <= min(len1 + max_indels, len2); ++i2) {
        auto ind = i2 + max_indels - len1;
        base[ind] = lizard_tail(static_cast<int>(len2 - i2));
    }

    vector<int> new_base(2*max_indels + 1, -INF);
    for (int i1 = static_cast<int>(len1) - 1; i1 >= 0; --i1) {
        for (int i2 = i1 + max_indels; i2 >= i1 - max_indels; --i2) {
            int inx = max_indels + i2 - i1;
            if ((i2 < 0) || (i2 > len2)) {
                new_base[inx] = -INF;
            } else if (i2 == len2) {
                new_base[inx] = lizard_tail(static_cast<int>(len1) - i1);
            } else {
                int m = (s1[i1] == s2[i2]) ? match : mismatch;

                if (inx == 0)
                    new_base[inx] = max( { m + base[inx], indel + new_base[inx + 1]} );
                else if (inx == 2*max_indels)
                    new_base[inx] = max( { m + base[inx], indel + base[inx - 1] } );
                else
                    new_base[inx] = max( { m + base[inx], indel + base[inx - 1], indel + new_base[inx + 1]} );
            }
        }

        std::swap(base, new_base);
    }

    return max(-INF, base[max_indels]); // Score is always finite due to possibility to align using mismatches
}

// vim: ts=4:sw=4
