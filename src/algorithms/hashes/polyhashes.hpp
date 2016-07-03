#pragma once

#include <vector>
#include <seqan/sequence.h>

namespace algorithms {
    // TODO Reimplement it as "generator" with O(1) memory consumption
    // WARNING: require seqan::length() for string s
    // Make sure that function seqan::string(T s) is defined
    template<typename T>
    std::vector<size_t> polyhashes(const T &s, size_t K) {
        if (seqan::length(s) < K) {
            return {  };
        }
        const size_t p = 7;
        std::vector<size_t> result(seqan::length(s) - K + 1);

        size_t first = 0;
        size_t p_pow_K = 1;
        for (size_t i = 0; i < K; ++i) {
            first *= p;
            first += unsigned(s[i]);
            p_pow_K *= p;
        }

        result[0] = first;
        for (size_t i = 1; i < result.size(); ++i) {
            first *= p;
            first += unsigned(s[K + i - 1]); //size_t(s[K + i - 1]);
            first -= unsigned(s[i - 1]) * p_pow_K; //size_t(s[i - 1]) * p_pow_K;
            result[i] = first;
        }
        return result;
    }
}