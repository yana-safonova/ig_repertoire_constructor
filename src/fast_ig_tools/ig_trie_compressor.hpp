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

template <typename Tletter = seqan::Dna5>
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

    template <typename Tcont>
    Trie(const Tcont &cont) : Trie() {
        for (const auto &s : cont) {
            add(s);
        }
    }

    template <typename T>
    void add(const T &s) {
        assert(!isCompressed());

        TrieNode *p = root;

        for (size_t i = 0; i < seqan::length(s); ++i) {
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
            root->compress_to_shortest();
        }
    }

    std::vector<size_t> checkout_vector() {
        if (!isCompressed()) {
            compress();
        }

        std::vector<size_t> result(this->size());
        root->checkout_vector(result);

        return result;

    }

    std::unordered_map<size_t, size_t> checkout(size_t nbucket = 0) {
        if (nbucket == 0) {
            nbucket = root->leaves_count();
        }

        if (!isCompressed()) {
            compress();
        }

        std::unordered_map<size_t, size_t> result(nbucket);
        root->checkout(result);

        return result;
    }

    std::unordered_map<size_t, std::vector<size_t>> checkout_ids(size_t nbucket = 0) {
        if (nbucket == 0) {
            nbucket = root->leaves_count();
        }

        if (!isCompressed()) {
            compress();
        }

        std::unordered_map<size_t, std::vector<size_t>> result(nbucket);
        root->checkout(result);

        return result;
    }

private:
    class TrieNode;
    static constexpr size_t card = seqan::ValueSize<Tletter>::VALUE;

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
        void compress_to_shortest(TrieNode *p = nullptr) {
            if (!p && ids) {
                p = this;
            }

            if (ids) {
                ids->target_node = p;
            }

            for (auto &child : children) {
                if (child) {
                    child->compress_to_shortest(p);
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


        void checkout(std::unordered_map<size_t, size_t> &result) const {
            if (ids && !ids->empty()) {
                size_t id = ids->target_node->ids->represent();
                result[id] += ids->size();
            }

            // DFS
            for (const auto &child : children) {
                if (child) {
                    child->checkout(result);
                }
            }
        }

        void checkout(std::unordered_map<size_t, std::vector<size_t>> &result) const {
            if (ids && !ids->empty()) {
                assert(!ids->target_node->ids->empty());
                size_t id = ids->target_node->ids->represent();
                result[id].insert(result[id].end(), ids->id_vector.cbegin(), ids->id_vector.cend());
            }

            // DFS
            for (const auto &child : children) {
                if (child) {
                    child->checkout(result);
                }
            }
        }

        size_t leaves_count() const {
            size_t result = 0;
            bool is_leaf = true;

            for (const auto &child : children) {
                if (child) {
                    is_leaf = false;
                    result += child->leaves_count();
                }
            }

            if (is_leaf) {
                result += 1;
            }

            return result;
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

template<typename Tletter=seqan::Dna5>
std::vector<size_t> compressed_reads_indices(const std::vector<seqan::String<Tletter>> &reads) {
    Trie<seqan::Dna5> trie(reads);

    return trie.checkout_vector();
}

template<typename Tletter=seqan::Dna5>
std::vector<seqan::String<Tletter>> compressed_reads(const std::vector<seqan::String<Tletter>> &reads) {
    Trie<seqan::Dna5> trie(reads);

    auto indices = trie.checkout_vector();

    std::vector<seqan::String<Tletter>> result;
    for (size_t i = 0; i < indices.size(); ++i) {
        if (indices[i] == i) {
            result.push_back(reads[i]);
        }
    }
}
}
// vim: ts=4:sw=4
