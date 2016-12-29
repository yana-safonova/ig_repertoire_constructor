#pragma once

#include <cassert>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <exception>
#include <verify.hpp>

#include <seqan/seq_io.h>
#include "fast_ig_tools.hpp"
using seqan::length;


// TODO Reimplement it as "generator" with O(1) memory consumption
template<typename T>
std::vector<size_t> polyhashes(const T &s, size_t K) {
    if (length(s) < K) {
        return {  };
    }

    const size_t p = 7;
    std::vector<size_t> result(length(s) - K + 1);

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

template<typename T1, typename T2 = T1>
int hamming_rtrim(const T1& s1, const T2 &s2) {
    size_t len = std::min<size_t>(length(s1), length(s2));

    int res = 0;
    for (size_t i = 0; i < len; ++i) {
        res += (s1[i] != s2[i]) ? 1 : 0;
    }

    return res;
}


template<typename T>
void remove_duplicates(std::vector<T> &v, bool sorted = false) {
    // Remove duplicated items
    if (!sorted) {
        std::sort(v.begin(), v.end());
    }
    auto it = std::unique_copy(v.cbegin(), v.cend(), v.begin());
    v.resize(std::distance(v.begin(), it));
}


using KmerIndex = std::unordered_map<size_t, std::vector<size_t>>;


template<typename T>
KmerIndex kmerIndexConstruction(const std::vector<T> &input_reads, size_t K) {
    size_t initial_hashtable_size = (K <= 13) ? (1 << (2*K)) : (input_reads.size() * 200);
    KmerIndex kmer2reads(initial_hashtable_size);

    for (size_t j = 0; j < input_reads.size(); ++j) {
        for (size_t hash : polyhashes(input_reads[j], K)) {
            kmer2reads[hash].push_back(j); // Already sorted. Nice!
        }
    }

    for (auto &kv : kmer2reads) {
        remove_duplicates(kv.second, true);
    }

    return kmer2reads;
}


template<typename T>
size_t count_unique(std::vector<T> v) {
   remove_duplicates(v);

   return v.size();
}


template<typename T>
std::vector<size_t> find_candidates(const T &read,
                                    const KmerIndex &kmer2reads,
                                    size_t target_size,
                                    unsigned tau, size_t K,
                                    unsigned strategy) {
    size_t required_read_length = (strategy != 0) ? (K * (tau + strategy)) : 0;
    if (length(read) < required_read_length) {
        return {  };
    }

    std::vector<size_t> cand;

    if (strategy == 0) { // Simple O(N*M) strategy
        cand.resize(target_size);
        std::iota(cand.begin(), cand.end(), 0);
    } else { // Minimizers strategy
        std::vector<size_t> multiplicities;

        auto hashes = polyhashes(read, K);

        multiplicities.reserve(hashes.size());
        for (size_t hash : hashes) {
            auto it = kmer2reads.find(hash);
            if (it != kmer2reads.cend()) {
                multiplicities.push_back(it->second.size());
            } else {
                multiplicities.push_back(0);
            }
        }

        std::vector<size_t> ind = optimal_coverage(multiplicities, K, tau + strategy);

        std::unordered_map<size_t, size_t> hits;
        for (size_t i : ind) {
            size_t hash = hashes[i];
            auto it = kmer2reads.find(hash);
            if (it != kmer2reads.cend()) {
                for (size_t candidate_index : it->second) {
                    ++hits[candidate_index];
                }
            }
        }

        for (const auto &kv : hits) {
            if (kv.second >= strategy) {
                cand.push_back(kv.first);
            }
        }
    }

    return cand;
}


template<typename T, typename Tf>
Graph tauDistGraph(const std::vector<T> &input_reads,
                   const KmerIndex &kmer2reads,
                   const Tf &dist_fun,
                   unsigned tau,
                   unsigned K,
                   unsigned strategy,
                   size_t &num_of_dist_computations) {
    Graph g(input_reads.size());

    std::atomic<size_t> atomic_num_of_dist_computations;
    atomic_num_of_dist_computations = 0;

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 8))
    for (size_t j = 0; j < input_reads.size(); ++j) {
        auto cand = find_candidates(input_reads[j], kmer2reads, input_reads.size(), tau, K, strategy);

        size_t len_j = length(input_reads[j]);

        for (size_t i : cand) {
            size_t len_i = length(input_reads[i]);
            if (len_j < len_i || (len_i == len_j && j < i)) {
                size_t dist = dist_fun(input_reads[j], input_reads[i]);

                atomic_num_of_dist_computations += 1;

                if (dist <= tau) {
                    g[j].push_back( { i, dist } );
                }
            }
        }
    }

    // Undirecting
    auto gg = g;
    for (size_t i = 0; i < gg.size(); ++i) {
        for (const auto &_ : gg[i]) {
            g[_.first].push_back( { i, _.second } );
        }
    }
    gg.clear(); // Free memory

    SEQAN_OMP_PRAGMA(parallel for schedule(guided, 8))
    for (size_t j = 0; j < g.size(); ++j) {
        remove_duplicates(g[j]);
    }

    num_of_dist_computations = atomic_num_of_dist_computations;

    return g;
}


template<typename T, typename Tf>
Graph tauMatchGraph(const std::vector<T> &input_reads,
                    const std::vector<T> &reference_reads,
                    const KmerIndex &kmer2reads,
                    const Tf &dist_fun,
                    unsigned tau,
                    unsigned K,
                    unsigned strategy,
                    size_t &num_of_dist_computations) {
    Graph g(input_reads.size());

    std::atomic<size_t> atomic_num_of_dist_computations;
    atomic_num_of_dist_computations = 0;

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 8))
    for (size_t j = 0; j < input_reads.size(); ++j) {
        auto cand = find_candidates(input_reads[j], kmer2reads, reference_reads.size(), tau, K, strategy);

        for (size_t i : cand) {
            unsigned dist = dist_fun(input_reads[j], reference_reads[i]);

            atomic_num_of_dist_computations += 1;

            if (dist <= tau) {
                g[j].push_back( { i, dist } );
            }
        }
    }

    num_of_dist_computations = atomic_num_of_dist_computations;

    return g;
}

// vim: ts=4:sw=4
