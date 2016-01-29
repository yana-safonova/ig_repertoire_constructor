#pragma once

#include <cassert>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <exception>

#include <seqan/seq_io.h>
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
        first += size_t(s[i]);
        p_pow_K *= p;
    }

    result[0] = first;
    for (size_t i = 1; i < result.size(); ++i) {
        first *= p;
        first += size_t(s[K + i - 1]);
        first -= size_t(s[i - 1]) * p_pow_K;
        result[i] = first;
    }

    return result;
}


// TODO cover by tests
// and then, refactor it!!!!!!!!111111111111oneoneone
std::vector<size_t> optimal_coverage(const std::vector<int> &costs, size_t K, size_t ss_len = 3) {
    assert(ss_len >= 1);
    assert(costs.size() + K - 1 >= ss_len * K);

    const int INF = 1 << 30; // TODO Use exact value

    std::vector<std::vector<int>> icosts(ss_len, std::vector<int>(costs.size()));

    // Fill by cummin
    icosts[0][0] = costs[0];
    for (size_t i = 1; i < icosts[0].size(); ++i) {
        icosts[0][i] = std::min(costs[i], icosts[0][i - 1]);
    }

    for (size_t j = 1; j < ss_len; ++j) { // ss_len == 1 is useless
        // Kill first K*j elements
        for (size_t i = 0; i < K*j; ++i) {
            icosts[j][i] = INF;
        }

        for (size_t i = K*j; i < icosts[j].size(); ++i) {
            icosts[j][i] = std::min(icosts[j][i - 1],
                                    costs[i] + icosts[j - 1][i - K]);
        }
    }

    auto ans = icosts[ss_len - 1][icosts[ss_len - 1].size() - 1];

    assert(ans < INF);

    std::vector<size_t> result(ss_len);
    // Backward reconstruction
    size_t i = icosts[ss_len - 1].size() - 1;
    size_t j = ss_len - 1;

    while (j > 0) {
        if (icosts[j][i] == costs[i] + icosts[j - 1][i - K]) { // Take i-th element
            result[j] = i;
            i -= K;
            j -= 1;
        } else {
            i -= 1;
        }
    }
    assert(j == 0);
    // Find first element
    size_t ii = i;
    while (icosts[0][i] != costs[ii]) {
        --ii;
    }
    result[0] = ii;


    // Checking
    assert(result.size() == ss_len);
    int sum = 0;
    for (size_t i : result) {
        sum += costs.at(i);
    }
    assert(ans == sum);

    return result;
}


template<typename T1, typename T2 = T1>
int hamming_rtrim(const T1& s1, const T2 &s2) {
    size_t len = std::min<size_t>(length(s1), length(s2));

    size_t res = 0;
    for (size_t i = 0; i < len; ++i) {
        res += (s1[i] != s2[i]) ? 1 : 0;
    }

    return res;
}


using KmerIndex = std::unordered_map<size_t, std::vector<size_t>>;


template<typename T>
KmerIndex kmerIndexConstruction(const std::vector<T> &input_reads, size_t K) {
    KmerIndex kmer2reads(input_reads.size() * 200); // Count them more careful

    for (size_t j = 0; j < input_reads.size(); ++j) {
        for (size_t hash : polyhashes(input_reads[j], K)) {
            kmer2reads[hash].push_back(j); // Already sorted. Nice!
        }
    }

    return kmer2reads;
}


enum Strategy { NAIVE, SINGLE, PAIR, TRIPLE };


const char* toCString(Strategy strategy) {
    switch (strategy) {
        case Strategy::NAIVE:
            return "Naive strategy";
        case Strategy::SINGLE:
            return "Basic <<single>> strategy";
        case Strategy::PAIR:
            return "<<Pair>> strategy";
        case Strategy::TRIPLE:
            return "<<Triple>> strategy";
        default:
            return "Unknown strategy";
    }
}


template<typename T>
void remove_duplicates(std::vector<T> &v) {
    // Remove duplicated items
    std::sort(v.begin(), v.end());
    auto it = std::unique_copy(v.cbegin(), v.cend(), v.begin());
    v.resize(std::distance(v.begin(), it));
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
                                    int tau, size_t K,
                                    Strategy strategy) {
    int strategy_int = int(strategy);
    size_t required_read_length = (strategy_int != 0) ? (K * (tau + strategy_int)) : 0;
    if (length(read) < required_read_length) {
        return {  };
    }

    std::vector<int> costs;

    auto hashes = polyhashes(read, K);

    costs.reserve(hashes.size());
    for (size_t hash : hashes) {
        auto it = kmer2reads.find(hash);
        if (it != kmer2reads.cend()) // TODO check it
            costs.push_back(it->second.size());
        else
            costs.push_back(0);
    }

    std::vector<size_t> cand;

    if (strategy == Strategy::NAIVE) { // Simple O(N*M) strategy
        cand.resize(target_size);
        std::iota(cand.begin(), cand.end(), 0);
    } else if (strategy == Strategy::SINGLE) { // basic <<single>> strategy
        std::vector<size_t> ind = optimal_coverage(costs, K, tau + 1);

        for (size_t i : ind) {
            size_t hash = hashes[i];
            auto it = kmer2reads.find(hash);
            if (it != kmer2reads.cend()) {
                std::copy(it->second.cbegin(), it->second.cend(), std::back_inserter(cand));
            }
        }
    } else if (strategy == Strategy::PAIR) { // <<pair>> strategy
        std::vector<size_t> ind = optimal_coverage(costs, K, tau + 2);

        std::vector<const std::vector<size_t>*> reads_set;

        reads_set.reserve(ind.size());
        for (size_t i : ind) {
            size_t hash = hashes[i];
            auto it = kmer2reads.find(hash);
            if (it != kmer2reads.cend()) {
                reads_set.push_back(&(it->second));
            }
        }

        for (size_t i = 0; i < reads_set.size(); ++i) {
            for (size_t j = i + 1; j < reads_set.size(); ++j) {
                std::set_intersection(reads_set[i]->cbegin(), reads_set[i]->cend(),
                                      reads_set[j]->cbegin(), reads_set[j]->cend(),
                                      std::back_inserter(cand));
            }
        }
    } else if (strategy == Strategy::TRIPLE) { // <<triple>> strategy (fastest, but with complicated preprocessing)
        std::vector<size_t> ind = optimal_coverage(costs, K, tau + 3);

        std::vector<const std::vector<size_t>*> reads_set;

        reads_set.reserve(ind.size());
        for (size_t i : ind) {
            size_t hash = hashes[i];
            auto it = kmer2reads.find(hash);
            if (it != kmer2reads.cend()) {
                reads_set.push_back(&(it->second));
            }
        }

        for (size_t i = 0; i < reads_set.size(); ++i) {
            for (size_t j = i + 1; j < reads_set.size(); ++j) {
                for (size_t k = j + 1; k < reads_set.size(); ++k) {
                    std::vector<size_t> tmp;
                    tmp.reserve(reads_set[i]->size() + reads_set[j]->size());

                    std::set_intersection(reads_set[i]->cbegin(), reads_set[i]->cend(),
                                          reads_set[j]->cbegin(), reads_set[j]->cend(),
                                          std::back_inserter(tmp));
                    std::set_intersection(tmp.cbegin(), tmp.cend(),
                                          reads_set[k]->cbegin(), reads_set[k]->cend(),
                                          std::back_inserter(cand));
                }
            }
        }
    } else {
        throw std::invalid_argument("Unknown strategy");
    }

    // Remove duplicated items
    remove_duplicates(cand);

    return cand;
}


using Graph = std::vector<std::vector<std::pair<size_t, int>>>;


size_t numEdges(const Graph &graph) {
    size_t nE = 0;

    for (size_t i = 0; i < graph.size(); ++i) {
        for (const auto &neib : graph[i]) {
            if (i != neib.first) nE += 1; // Exclude loop edges
        }
    }
    nE /= 2;

    return nE;
}


void write_metis_graph(const Graph &graph,
                       const std::string &filename) {
    std::ofstream out(filename);

    // Count the numder of vertices and the number of edges
    size_t nV = graph.size();
    size_t nE = numEdges(graph);

    out << nV << " " << nE << " 001\n"; // See http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf

    for (size_t i = 0; i < graph.size(); ++i) {
        for (const auto &neib : graph[i]) {
            if (i == neib.first) continue; // Exclude loop edges
            out << neib.first + 1 << " " << neib.second << " ";
        }
        out << "\n";
    }
}


void write_metis_graph(const Graph &graph,
                       const std::vector<size_t> &weights,
                       const std::string &filename) {
    std::ofstream out(filename);

    // Count the numder of vertices and the number of edges
    size_t nV = graph.size();
    size_t nE = numEdges(graph);

    out << nV << " " << nE << " 011\n"; // See http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf

    for (size_t i = 0; i < graph.size(); ++i) {
        out << weights[i] << " ";
        for (const auto &neib : graph[i]) {
            if (i == neib.first) continue; // Exclude loop edges
            out << neib.first + 1 << " " << neib.second << " ";
        }
        out << "\n";
    }
}

// vim: ts=4:sw=4
