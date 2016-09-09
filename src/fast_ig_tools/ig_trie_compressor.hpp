#pragma once

#include <algorithm>
#include <array>
#include <boost/pool/object_pool.hpp>
#include <cassert>
#include <limits>
#include <memory>
#include <unordered_map>
#include <vector>

#include <seqan/seq_io.h>

namespace fast_ig_tools {

template <typename TValue = seqan::Dna5>
class Trie {
public:
    Trie() {
        // Do not use :member initialization syntax here. pool should be initialized BEFORE root
        root = pool.construct();
    }
    Trie(const Trie &) = delete;
    Trie &operator=(const Trie &) = delete;
    Trie(Trie &&) = default;
    Trie &operator=(Trie &&) = default;
    ~Trie() = default;


    size_t size() const {
        return size__;
    }

    template <typename TCont>
    Trie(const TCont &cont) : Trie(cont.cbegin(), cont.cend()) { }

    template <typename TIter>
    Trie(TIter b, TIter e) : Trie() {
        for (; b !=e; ++b) {
            add(*b);
        }
    }

    template <typename T>
    void add(const T &s) {
        assert(!isCompressed());

        TrieNode *p = root;

        for (size_t i = 0; i < seqan::length(s); ++i) {
			if (p->ids) {
				// Ok, nice, we have found a sequence that is a prefix of the current read
				// Join current sequence to the found prefix
				// Do not construct trie further
				break;
			}
            size_t el = seqan::ordValue(s[i]);
            assert((0 <= el) && (el < p->children.size()));

            if (!p->children[el]) {
                p->children[el] = pool.construct();
            }

            p = p->children[el];
        }

        if (!p->ids) {
            p->ids = new IdCounter;
        }

        p->ids->add(size__);
        size__++;
    }

    bool isCompressed() const {
        return compressed;
    }

    void compress() {
        if (!isCompressed()) {
            root->compress_to_prefix();
        }
    }

    std::vector<size_t> checkout() {
        if (!isCompressed()) {
            compress();
        }

        std::vector<size_t> result(this->size());
        root->checkout_vector(result);

        return result;

    }

private:
    class TrieNode;
    static constexpr size_t card = seqan::ValueSize<TValue>::VALUE;

    struct IdCounter {
        std::vector<size_t> id_vector;
        TrieNode *target_node = nullptr;

        void add(size_t id) {
            id_vector.push_back(id);
        }

        bool empty() const {
            return id_vector.empty();
        }

        size_t size() const {
            return id_vector.size();
        }

        size_t represent() const {
            assert(!empty());

            return id_vector[0];
        }
    };

    class TrieNode {
    public:
        static const size_t INF = std::numeric_limits<size_t>::max();
        std::array<TrieNode*, card> children;

        TrieNode() : ids{nullptr} {
            children.fill(nullptr);
        }

        IdCounter *ids;

        // Compression to prefix
        void compress_to_prefix(TrieNode *p = nullptr) {
            if (!p && ids) {
                p = this;
            }

            if (ids) {
                ids->target_node = p;
            }

            for (auto &child : children) {
                if (child) {
                    child->compress_to_prefix(p);
                }
            }
        }

        void checkout_vector(std::vector<size_t> &result) const {
            if (ids && !ids->empty()) {
                size_t target_id = ids->target_node->ids->represent();
                for (size_t id : ids->id_vector) {
                    result[id] = target_id;
                }
            }

            // DFS
            for (const auto &child : children) {
                if (child) {
                    child->checkout_vector(result);
                }
            }
        }

        ~TrieNode() {
            if (ids) {
                delete ids;
            }

            // Do not do this! Destructors will be called by object_pool
            // for (auto &child : children) {
            //     if (child)
            //         delete child;
            // }
        }
    };

    boost::object_pool<TrieNode> pool;
    TrieNode *root;
    bool compressed = false;
    size_t size__ = 0;
};

template <typename T>
using Decay = typename std::decay<T>::type;

template <typename TArray>
using ValueType = Decay<decltype((Decay<TArray>())[0])>;

template <typename TIter>
std::vector<size_t> compressed_reads_indices(TIter b, TIter e) {
    using TValue = ValueType<decltype(*(TIter()))>;
    Trie<TValue> trie(b, e);

    return trie.checkout();
}

template <typename TVector>
std::vector<size_t> compressed_reads_indices(const TVector &reads) {
    return compressed_reads_indices(reads.cbegin(), reads.cend());
}

template <typename TIter>
auto compressed_reads(TIter b, TIter e) -> std::vector<Decay<decltype(*(TIter()))>>  {
    auto indices = compressed_reads_indices(b, e);

    std::vector<Decay<decltype(*(TIter()))>> result;
    for (size_t i = 0; i < indices.size(); ++i, ++b) {
        if (indices[i] == i) {
            result.push_back(*b);
        }
    }

    return result;
}

template <typename TVector>
auto compressed_reads(const TVector &reads) -> std::vector<ValueType<TVector>>  {
    return compressed_reads(reads.cbegin(), reads.cend());
}
}
// vim: ts=4:sw=4
