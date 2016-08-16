#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <exception>
#include <fstream>
#include <limits>
#include <unordered_map>
#include <vector>

#include <verify.hpp>

#include <folly/SmallLocks.h>

#include "bitread.hpp"
#include <boost/unordered_map.hpp>
#include <seqan/seq_io.h>
#undef NDEBUG
#include <cassert>

namespace hogc {
template <typename TString>
std::vector<size_t> rolling_hashes(const TString &s, const size_t K, const size_t p = 7) {
    using seqan::length;
    using seqan::ordValue;

    if (length(s) < K) {
        return {};
    }

    std::vector<size_t> result(length(s) - K + 1);

    size_t first = 0;
    size_t p_pow_K = 1;
    for (size_t i = 0; i < K; ++i) {
        first *= p;
        first += ordValue(s[i]);
        p_pow_K *= p;
    }

    result[0] = first;
    for (size_t i = 1; i < result.size(); ++i) {
        first *= p;
        first += ordValue(s[K + i - 1]);        //size_t(s[K + i - 1]);
        first -= ordValue(s[i - 1]) * p_pow_K;  //size_t(s[i - 1]) * p_pow_K;
        result[i] = first;
    }

    return result;
}

template <typename TElement>
void remove_duplicates(std::vector<TElement> &v, bool sorted = false) {
    // Remove duplicated items
    if (!sorted) {
        std::sort(v.begin(), v.end());
    }
    auto it = std::unique_copy(v.cbegin(), v.cend(), v.begin());
    v.resize(std::distance(v.begin(), it));
}

template <typename TVertexWeight = size_t, typename TEdgeWeight = int>
class Graph {
public:
    Graph(size_t V, bool undirected = false) : edges__(V),
                                               undirected__{undirected},
                                               weighted_vertices__{false} {}
    Graph(const std::vector<TVertexWeight> &vw, bool undirected = false) : vertices__(vw),
                                                                           edges__(vw.size()),
                                                                           undirected__{undirected},
                                                                           weighted_vertices__{true} {}
    Graph(std::vector<TVertexWeight> &&vw, bool undirected = false) : vertices__(vw),
                                                                      edges__(vw.size()),
                                                                      undirected__{undirected},
                                                                      weighted_vertices__{true} {}

    void setVerticesWeights(std::vector<TVertexWeight> &&vw) {
        assert(edges__.size() == vw.size());
        vertices__ = vw;
        weighted_vertices__ = true;
    }

    void setVerticesWeights(const std::vector<TVertexWeight> &vw) {
        auto tmp = vw;
        setVerticesWeights(vw);
    }

    void normalize() {
        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 128))
        for (size_t i = 0; i < edges__.size(); ++i) {
            remove_duplicates(edges__[i]);
        }
    }

    bool isWeighedVertices() const {
        return weighted_vertices__;
    }

    const std::vector<std::vector<std::pair<size_t, TEdgeWeight>>> &edges() const {
        return edges__;
    }

    bool verify() const {
        if (isWeighedVertices()) {
            if (vertices__.size() != edges__.size()) {
                return false;
            }
        }

        // Check ordering & uniqueness
        for (const auto &neibs : edges__) {
            if (!std::is_sorted(neibs.cbegin(), neibs.cend())) {
                return false;
            }

            if (std::adjacent_find(neibs.cbegin(), neibs.cend()) != neibs.cend()) {
                return false;
            }
        }

        // Check symmetry
        if (undirected__) {
            for (size_t i = 0; i < edges__.size(); ++i) {
                for (const auto &neib : edges__[i]) {
                    size_t j = neib.first;

                    if (i == j) {
                        return false;
                    }

                    auto inv_edge = std::make_pair(i, neib.second);
                    if (!std::binary_search(edges__[j].cbegin(), edges__[j].cend(), inv_edge)) {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    size_t nVertices() const {
        return edges__.size();
    }

    size_t nEdges() const {
        size_t nE = 0;

        for (const auto &neibs : edges__) {
            nE += neibs.size();
        }

        if (isUndirected()) {
            nE /= 2;
        }

        return nE;
    }

    bool isUndirected() const {
        return undirected__;
    }

    // TODO Use merge inst of sorting (set_union)
    void undirect() {
        assert(!isUndirected());
        auto edges_copy = edges__;
        for (size_t i = 0; i < edges_copy.size(); ++i) {
            for (const auto &edge : edges_copy[i]) {
                edges__[edge.first].push_back({i, edge.second});
            }
        }
        edges_copy.clear();  // Free memory

        SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 128))
        for (size_t j = 0; j < edges__.size(); ++j) {
            remove_duplicates(edges__[j]);
        }

        undirected__ = true;
    }

    void writeMetis(const std::string &filename) const {
        assert(isUndirected());

        std::ofstream out(filename);

        if (isWeighedVertices()) {
            out << nVertices() << " " << nEdges() << " 011\n";  // See http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf

            for (size_t i = 0; i < edges__.size(); ++i) {
                out << vertices__[i] << " ";
                for (const auto &edge : edges__[i]) {
                    out << edge.first + 1 << " " << edge.second << " ";
                }
                out << "\n";
            }
        } else {
            out << nVertices() << " " << nEdges() << " 001\n";  // See http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf

            for (const auto &edges : edges__) {
                for (const auto &edge : edges) {
                    out << edge.first + 1 << " " << edge.second << " ";
                }
                out << "\n";
            }
        }
    }

    void addEdge(size_t i, size_t j, TEdgeWeight w) {
        // Thread-safe iff i1 != i2
        edges__[i].push_back({j, w});
    }

private:
    std::vector<TVertexWeight> vertices__;
    std::vector<std::vector<std::pair<size_t, TEdgeWeight>>> edges__;
    bool undirected__;
    bool weighted_vertices__;
};

bool check_repr_kmers_consistancy(const std::vector<size_t> &answer,
                                  const std::vector<int> &multiplicities,
                                  size_t K, size_t n) {
    if (!std::is_sorted(answer.cbegin(), answer.cend())) {
        return false;
    }

    // Check answer size
    if (answer.size() != n) {
        return false;
    }

    // Check k-mer overlapping
    for (size_t i = 1; i < answer.size(); ++i) {
        if (answer[i] - answer[i - 1] < K) {
            return false;
        }
    }

    // K-mers should belong interval
    for (size_t kmer_i : answer) {
        if (kmer_i >= multiplicities.size()) {
            return false;
        }
    }

    return true;
}

// TODO cover by tests
// and then, refactor it!!!!!!!!111111111111oneoneone
// TODO Rename it to be consistent with the paper
// TODO int -> unsigned or size_t
std::vector<size_t> optimal_coverage(const std::vector<int> &multiplicities,
                                     size_t K, size_t n = 3) {
    assert(n >= 1);
    assert(multiplicities.size() + K - 1 >= n * K);

    const int INF = 1 << 30;  // TODO Use exact value

    std::vector<std::vector<int>> imults(n, std::vector<int>(multiplicities.size()));

    // Fill by cummin
    imults[0][0] = multiplicities[0];
    for (size_t i = 1; i < imults[0].size(); ++i) {
        imults[0][i] = std::min(multiplicities[i], imults[0][i - 1]);
    }

    for (size_t j = 1; j < n; ++j) {  // n == 1 is useless
        // Kill first K*j elements
        for (size_t i = 0; i < K * j; ++i) {
            imults[j][i] = INF;
        }

        for (size_t i = K * j; i < imults[j].size(); ++i) {
            imults[j][i] = std::min(imults[j][i - 1],
                                    multiplicities[i] + imults[j - 1][i - K]);
        }
    }

    auto ans = imults[n - 1][imults[n - 1].size() - 1];

    VERIFY(ans < INF);

    std::vector<size_t> result(n);
    // Backward reconstruction
    size_t i = imults[n - 1].size() - 1;
    size_t j = n - 1;

    while (j > 0) {
        if (imults[j][i] == multiplicities[i] + imults[j - 1][i - K]) {  // Take i-th element
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
    while (imults[0][i] != multiplicities[ii]) {
        --ii;
    }
    result[0] = ii;

    // Checking
    int sum = 0;
    for (size_t i : result) {
        sum += multiplicities.at(i);
    }
    assert(ans == sum);

    assert(check_repr_kmers_consistancy(result, multiplicities, K, n));

    return result;
}

template <typename T>
class locked_vector : public std::vector<T> {
public:
    //  Inherit constructors
    using std::vector<T>::vector;

    //Thread-safe push_back()
    template <typename... Args>
    void push_back(Args &&... args) {
        folly::MSLGuard guard(lock__);
        std::vector<T>::push_back(std::forward<Args>(args)...);
    }

private:
    folly::MicroSpinLock lock__;
};

// Inspired by folly::SpinLockArray, but not POD
template <class T, size_t CACHE_LINE_SIZE = 64>
class SpinLockVector final {
public:
    SpinLockVector(size_t N) : data_(N) {}
    void resize(size_t N) {
        data_.resize(N);
    }

    T &operator[](size_t i) {
        return data_[i].lock;
    }

    const T &operator[](size_t i) const {
        return data_[i].lock;
    }

    size_t size() const { return data_.size(); }

private:
    struct PaddedSpinLock {
        PaddedSpinLock() : lock() {}
        T lock;
        char padding[CACHE_LINE_SIZE - sizeof(T)];
    };
    static_assert(sizeof(PaddedSpinLock) == CACHE_LINE_SIZE,
                  "Invalid size of PaddedSpinLock");

    // Check if T can theoretically cross a cache line.
    // NOTE: It should be alignof(std::max_align_t), but max_align_t
    // isn't supported by gcc 4.6.2.
    static_assert(alignof(MaxAlign) > 0 &&
                          CACHE_LINE_SIZE % alignof(MaxAlign) == 0 &&
                          sizeof(T) <= alignof(MaxAlign),
                  "T can cross cache line boundaries");

    char padding_[CACHE_LINE_SIZE];
    std::vector<PaddedSpinLock> data_;
} __attribute__((aligned));

template <typename TIIter, typename TOIter>
TOIter freq_filter(TIIter b, TIIter e, TOIter out, size_t limit = 1) {
    if (b == e) {
        return out;
    }

    TIIter cur = b, prev = b;
    ++cur;
    size_t cur_freq = 1;

    while (cur != e) {
        assert(*prev <= *cur);

        if (*cur == *prev) {
            ++cur_freq;
        } else {
            if (cur_freq >= limit) {
                *out++ = *prev;
            }
            cur_freq = 1;
        }

        ++prev;
        ++cur;
    }

    if (cur_freq >= limit) {
        *out++ = *prev;
    }

    return out;
}

template <typename TIndex = uint32_t, typename TPos = uint32_t>
class KmerIndex {
public:
    template <typename TString>
    KmerIndex(const std::vector<TString> &reads,
              size_t K,
              size_t bins = 1 << 20,
              size_t p = 7) : K{K}, p{p}, N{reads.size()}, bins{bins}, positional_index__(bins) {
        assert(std::numeric_limits<TIndex>::max() >= N);

        for (size_t i = 0; i < reads.size(); ++i) {
            lengths__.push_back(seqan::length(reads[i]));
        }

        SpinLockVector<folly::MicroSpinLock> locks(bins);
        positional_index__.resize(bins);

        SEQAN_OMP_PRAGMA(parallel for)
        for (size_t i = 0; i < reads.size(); ++i) {
            auto hashes = rolling_hashes__(reads[i]);

            for (size_t pos = 0; pos < hashes.size(); ++pos) {
                size_t bin = hashes[pos];
                folly::MSLGuard guard(locks[bin]);
                positional_index__[bin][pos].push_back(i);
            }
        }

        SEQAN_OMP_PRAGMA(parallel for)
        for (size_t bin = 0; bin < bins; ++bin) {
            for (auto &kv : positional_index__[bin]) {
                auto &indices = kv.second;
                indices.shrink_to_fit();
                // std::sort(indices.begin(), indices.end());
            }
        }
    }

    template <typename TString>
    std::vector<size_t> candidates(const TString &read,
                                   size_t tau,
                                   size_t strategy = 1,
                                   size_t shift = 0,
                                   bool len_filter = true) const {
        size_t required_read_length = (strategy != 0) ? (K * (tau + strategy)) : 0;
        size_t len = seqan::length(read);
        if (len < required_read_length + shift) {
            return {};
        }

        std::vector<size_t> cand;

        if (strategy == 0) {  // Simple O(N*M) strategy
            cand.resize(N);
            std::iota(cand.begin(), cand.end(), 0);
        } else {  // Minimizers strategy
            auto hashes = rolling_hashes__(read);

            // TODO Use size_t inst of int
            std::vector<int> multiplicities;
            std::vector<const std::vector<TIndex> *> responses;
            for (size_t i = 0; i + shift < hashes.size(); ++i) {
                auto response = &query__(hashes[i + shift], i);
                multiplicities.push_back(response->size());
                responses.push_back(response);
            }

            std::vector<size_t> ind = optimal_coverage(multiplicities, K, tau + strategy);

            for (TPos i : ind) {
                for (TIndex candidate_index : *responses[i]) {
                    cand.push_back(candidate_index);
                }
            }

            std::sort(cand.begin(), cand.end());
            auto it = freq_filter(cand.cbegin(), cand.cend(), cand.begin(), strategy);
            cand.resize(it - cand.begin());
        }

        if (len_filter) {
            auto it = std::copy_if(cand.cbegin(), cand.cend(), cand.begin(),
                                   [this, shift, len](size_t i) -> bool { return lengths__[i] + shift >= len; });
            cand.resize(it - cand.begin());
        }

        return cand;
    }

    template <typename TString>
    std::vector<size_t> candidates_negshift(const TString &read,
                                            size_t tau,
                                            size_t strategy = 1,
                                            size_t shift = 0,
                                            bool len_filter = true) const {
        size_t required_read_length = (strategy != 0) ? (K * (tau + strategy)) : 0;
        size_t len = seqan::length(read);
        if (len < required_read_length) {
            return {};
        }

        std::vector<size_t> cand;

        if (strategy == 0) {  // Simple O(N*M) strategy
            cand.resize(N);
            std::iota(cand.begin(), cand.end(), 0);
        } else {  // Minimizers strategy
            auto hashes = rolling_hashes__(read);

            // TODO Use size_t inst of int
            std::vector<int> multiplicities;
            std::vector<const std::vector<TIndex> *> responses;
            for (size_t i = 0; i < hashes.size(); ++i) {
                auto response = &query__(hashes[i], i + shift);
                multiplicities.push_back(response->size());
                responses.push_back(response);
            }

            std::vector<size_t> ind = optimal_coverage(multiplicities, K, tau + strategy);

            for (TPos i : ind) {
                for (TIndex candidate_index : *responses[i]) {
                    cand.push_back(candidate_index);
                }
            }

            std::sort(cand.begin(), cand.end());
            auto it = freq_filter(cand.cbegin(), cand.cend(), cand.begin(), strategy);
            cand.resize(it - cand.begin());
        }

        if (len_filter) {
            auto it = std::copy_if(cand.cbegin(), cand.cend(), cand.begin(),
                                   [this, shift, len](size_t i) -> bool { return lengths__[i] >= len + shift; });
            cand.resize(it - cand.begin());
        }

        return cand;
    }

    const size_t K, p, N, bins;

private:
    template <typename TString>
    std::vector<size_t> rolling_hashes__(const TString &read) const {
        std::vector<size_t> hashes = rolling_hashes(read, K, p);
        for (size_t &hash : hashes) {
            hash %= bins;
        }

        return hashes;
    }

    const std::vector<TIndex> &query__(size_t hash, TPos pos) const {
        static std::vector<TIndex> empty = {};
        auto it = positional_index__[hash].find(pos);
        return it != positional_index__[hash].cend() ? it->second : empty;
    }

    size_t count__(size_t hash, TPos pos) const {
        return query__(hash, pos).size();
    }

    std::vector<boost::unordered_map<TPos, std::vector<TIndex>>> positional_index__;
    std::vector<size_t> lengths__;
};

template <typename T1, typename T2 = T1>
int hamming_over(const T1 &s1, const T2 &s2, size_t shift = 0) {
    using seqan::length;
    if (length(s2) + shift < length(s1)) {
        return std::numeric_limits<int>::max() / 2;
    }

    size_t len = length(s1) - shift;

    int res = 0;
    for (size_t i = 0; i < len; ++i) {
        res += (s1[i + shift] != s2[i]) ? 1 : 0;
    }

    return res;
}

template <typename TString>
size_t discardedReads(const std::vector<TString> &reads,
                      size_t tau,
                      size_t K,
                      size_t strategy,
                      int overlap = 0) {
    if (overlap > 0) {
        if (strategy > 0 && (tau + strategy) * K > static_cast<size_t>(overlap)) {
            return reads.size();
        } else {
            return 0;
        }
    }

    size_t shift = static_cast<size_t>(-overlap);
    size_t min_read_len = (strategy > 0) ? shift + (tau + strategy) * K : shift;

    size_t count = 0;
    for (const auto &read : reads) {
        if (seqan::length(read) < min_read_len) {
            ++count;
        }
    }

    return count;
}

template <typename TString>
Graph<> __hammingOverlapGraph_old(const std::vector<TString> &reads,
                                  size_t tau,
                                  int overlap = 0,
                                  size_t K = 10,
                                  size_t strategy = 1,
                                  size_t *num_of_dist_computations = nullptr) {
    Graph<> g(reads.size());

    KmerIndex<> kmer_index(reads, K);

    std::atomic<size_t> atomic_num_of_dist_computations;
    atomic_num_of_dist_computations = 0;

    SEQAN_OMP_PRAGMA(parallel for)
    for (size_t i = 0; i < reads.size(); ++i) {
        size_t len = seqan::length(reads[i]);

        assert(std::abs(overlap) <= len);

        size_t maximal_shift = (overlap <= 0) ? static_cast<size_t>(-overlap) : static_cast<size_t>(len - overlap);

        for (size_t shift = 0; shift <= maximal_shift; ++shift) {
            auto cand = kmer_index.candidates(reads[i], tau, strategy, shift);

            for (size_t j : cand) {
                if (i != j) {
                    size_t dist = hamming_over(reads[i], reads[j], shift);

                    // atomic_num_of_dist_computations += 1;

                    if (dist <= tau) {
                        g.addEdge(i, j, dist);
                    }
                }
            }
        }
    }

    if (num_of_dist_computations) {
        *num_of_dist_computations = atomic_num_of_dist_computations;
    }

    g.normalize();

    return g;
}

template <typename TString>
Graph<> hammingOverlapGraph(const std::vector<TString> &reads,
                            size_t tau,
                            int overlap = 0,
                            size_t K = 10,
                            size_t strategy = 1,
                            size_t *num_of_dist_computations = nullptr) {
    using fast_ig_tools::shifted_bitread;

    Graph<> g(reads.size()), gt(reads.size());

    std::vector<shifted_bitread<640>> bitreads(reads.size());

    INFO("Converting reads to bit representation...");
    SEQAN_OMP_PRAGMA(parallel for)
    for (size_t i = 0; i < reads.size(); ++i) {
        assert(length(reads[i]) <= 640);
        bitreads[i] = shifted_bitread<640>(reads[i]);
    }
    INFO("Converting reads done");

    INFO("k-mer index construction...");
    KmerIndex<> kmer_index(reads, K);
    INFO("k-mer index constructed");

    size_t n_of_dist_comps = 0;
    SEQAN_OMP_PRAGMA(parallel for reduction(+:n_of_dist_comps) schedule(dynamic, 128))
    for (size_t i = 0; i < reads.size(); ++i) {
        size_t len = seqan::length(reads[i]);

        if (std::abs(overlap) > len) {
            continue;
        }

        size_t maximal_shift = (overlap <= 0) ? static_cast<size_t>(-overlap) : static_cast<size_t>(len - overlap);

        for (size_t shift = 0; shift <= maximal_shift; ++shift) {
            auto cand = kmer_index.candidates(reads[i], tau, strategy, shift);

            for (size_t j : cand) {
                if (i != j) {
                    size_t dist = bitreads[i].dist(bitreads[j], shift);

                    ++n_of_dist_comps;

                    if (dist <= tau) {
                        g.addEdge(i, j, dist);
                        // g.addEdge(i, j, shift);
                    }
                }
            }
        }

        size_t maximal_neg_shift = (overlap <= 0) ? static_cast<size_t>(-overlap) : 0;
        for (size_t shift = 1; shift <= maximal_neg_shift; ++shift) {
            auto cand = kmer_index.candidates_negshift(reads[i], tau, strategy, shift);

            for (size_t j : cand) {
                if (i != j) {
                    size_t dist = bitreads[i].dist(bitreads[j], -static_cast<int>(shift));

                    ++n_of_dist_comps;

                    if (dist <= tau) {
                        gt.addEdge(i, j, dist);
                        // gt.addEdge(i, j, shift);
                    }
                }
            }
        }
    }

    if (num_of_dist_computations) {
        *num_of_dist_computations = n_of_dist_comps;
    }

    for (size_t i = 0; i < reads.size(); ++i) {
        for (const auto &neib : gt.edges()[i]) {
            g.addEdge(neib.first, i, neib.second);
        }
    }

    g.normalize();

    return g;
}

template <typename TString>
Graph<> __hammingGraphTail_old(const std::vector<TString> &reads,
                               size_t tau,
                               unsigned K = 10,
                               unsigned strategy = 1,
                               size_t *num_of_dist_computations = nullptr) {
    Graph<> g(reads.size());

    KmerIndex<> kmer_index(reads, K);

    std::atomic<size_t> atomic_num_of_dist_computations;
    atomic_num_of_dist_computations = 0;

    SEQAN_OMP_PRAGMA(parallel for)
    for (size_t i = 0; i < reads.size(); ++i) {
        size_t len = seqan::length(reads[i]);
        auto cand = kmer_index.candidates(reads[i], tau, strategy, 0);

        for (size_t j : cand) {
            size_t len_j = seqan::length(reads[j]);
            if (len < len_j || (len == len_j && i < j)) {
                unsigned dist = hamming_over(reads[i], reads[j], 0);

                atomic_num_of_dist_computations += 1;

                if (dist <= tau) {
                    g.addEdge(i, j, dist);
                }
            }
        }
    }

    if (num_of_dist_computations) {
        *num_of_dist_computations = atomic_num_of_dist_computations;
    }

    g.normalize();
    g.undirect();

    return g;
}

template <typename TString>
Graph<> hammingGraphIgnoreTails(const std::vector<TString> &reads,
                                size_t tau,
                                size_t K = 10,
                                size_t strategy = 1,
                                size_t *num_of_dist_computations = nullptr) {
    using fast_ig_tools::bitread;

    Graph<> g(reads.size());

    std::vector<bitread<640>> bitreads(reads.size());

    INFO("Converting reads to bit representation...");
    SEQAN_OMP_PRAGMA(parallel for)
    for (size_t i = 0; i < reads.size(); ++i) {
        assert(length(reads[i]) <= 640);
        bitreads[i] = bitread<640>(reads[i]);
    }
    INFO("Converting reads done");

    INFO("k-mer index construction...");
    KmerIndex<> kmer_index(reads, K);
    INFO("k-mer index constructed");

    size_t n_of_dist_comps = 0;

    INFO("Distances computation...");
    SEQAN_OMP_PRAGMA(parallel for reduction(+:n_of_dist_comps) schedule(dynamic, 128))
    for (size_t i = 0; i < reads.size(); ++i) {
        auto cand = kmer_index.candidates(reads[i], tau, strategy, 0, false);

        for (size_t j : cand) {
            if (bitreads[i].length() < bitreads[j].length() || (bitreads[i].length() == bitreads[j].length() && i < j)) {
                size_t dist = bitreads[i].dist(bitreads[j]);

                ++n_of_dist_comps;

                if (dist <= tau) {
                    g.addEdge(i, j, dist);
                }
            }
        }
    }
    INFO("Distances computation done");

    if (num_of_dist_computations) {
        *num_of_dist_computations = n_of_dist_comps;
    }

    INFO("Graph normalization...");
    g.normalize();
    g.undirect();
    INFO("Graph normalization done");

    return g;
}

template <typename TString>
Graph<> hammingGraphEqualLength(const std::vector<TString> &reads,
                                size_t tau,
                                size_t K = 10,
                                size_t strategy = 1,
                                size_t *num_of_dist_computations = nullptr) {
    using fast_ig_tools::bitread;
    using seqan::length;

    Graph<> g(reads.size());

    // Split reads
    boost::unordered_map<size_t, std::vector<TString>> bins;
    boost::unordered_map<size_t, std::vector<size_t>> bins2read;

    for (size_t i = 0; i < reads.size(); ++i) {
        size_t len = length(reads[i]);
        bins[len].push_back(reads[i]);
        bins2read[len].push_back(i);
    }

    size_t all_n_of_dist_comps = 0;
    for (const auto &kv : bins) {
        size_t len = kv.first;
        const auto &reads = kv.second;
        const auto &map_to_initial = bins2read[len];

        KmerIndex<> kmer_index(reads, K);

        std::vector<bitread<640>> bitreads(reads.size());

        SEQAN_OMP_PRAGMA(parallel for)
        for (size_t i = 0; i < reads.size(); ++i) {
            assert(length(reads[i]) <= 640);
            bitreads[i] = bitread<640>(reads[i]);
        }

        size_t n_of_dist_comps = 0;
        SEQAN_OMP_PRAGMA(parallel for reduction(+:n_of_dist_comps) schedule(dynamic, 128))
        for (size_t i = 0; i < reads.size(); ++i) {
            auto cand = kmer_index.candidates(reads[i], tau, strategy, 0, false);

            for (size_t j : cand) {
                if (i < j) {
                    size_t dist = bitreads[i].dist_equal_len(bitreads[j]);

                    ++n_of_dist_comps;

                    if (dist <= tau) {
                        g.addEdge(map_to_initial[i], map_to_initial[j], dist);
                    }
                }
            }
        }
        all_n_of_dist_comps += n_of_dist_comps;
    }

    if (num_of_dist_computations) {
        *num_of_dist_computations = all_n_of_dist_comps;
    }

    INFO("Graph normalization...");
    g.normalize();
    g.undirect();
    INFO("Graph normalization done");

    return g;
}

template <typename TString>
Graph<> hammingGraph(const std::vector<TString> &reads,
                     size_t tau,
                     bool ignore_tails,
                     size_t K = 10,
                     size_t strategy = 1,
                     size_t *num_of_dist_computations = nullptr) {
    return ignore_tails ? hammingGraphIgnoreTails(reads, tau, K, strategy, num_of_dist_computations) : hammingGraphEqualLength(reads, tau, K, strategy, num_of_dist_computations);
}

}  // namespace hogc
