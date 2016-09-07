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

    template <typename Tcont>
    Trie(const Tcont &cont) : Trie() {
        size_t i = 0;
        for (const auto &s : cont) {
            add(s, i);
            ++i;
        }
    }

    template <typename TcontRead, typename TcontAbun>
    Trie(const TcontRead &cont_read, const TcontAbun &cont_abun) : Trie() {
        size_t i = 0;
        // Zip!!! I badly need zip!

        auto it_read = cont_read.cbegin(), end_read = cont_read.cend();
        auto it_abun = cont_abun.cbegin(), end_abun = cont_abun.cend();
        for (; it_read != end_read && it_abun != end_abun; ++it_read, ++it_abun) {
            add(*it_read, i, *it_abun);
            ++i;
        }

        assert(it_read == end_read && it_abun == end_abun);  // contairners should have the same length
    }

    template <typename T, typename Tf>
    void add(const T &s, size_t id, const Tf &toIndex, size_t abundance = 1) {
        assert(!isCompressed());

        typename TrieNode::pointer_type p = root;

        for (size_t i = 0; i < seqan::length(s); ++i) {
            size_t el = toIndex(s[i]);
            assert((0 <= el) && (el < p->children.size()));

            if (!p->children[el]) {
                p->children[el] = pool.construct();
            }

            p = p->children[el];
        }

        if (!p->ids) {
            p->ids = new IdCounter;
        }

        p->ids->add(id, abundance);
    }

    template <typename T>
    void add(const T &s, size_t id) {
        auto to_size_t = [](const Tletter &letter) -> size_t { return seqan::ordValue(letter); };
        add(s, id, to_size_t);
    }

    bool isCompressed() const {
        return compressed;
    }

    void compress() {
        if (!isCompressed()) {
            root->compress_to_shortest();
        }
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
    static constexpr size_t card = seqan::ValueSize<Tletter>::VALUE;

    struct IdCounter {
        std::vector<size_t> id_vector;
        size_t count = 0;

        void add(size_t id, size_t abundance = 1) {
            id_vector.push_back(id);
            count += abundance;
        }

        bool empty() const {
            return count == 0;
        }

        size_t size() const {
            return count;
        }

        size_t represent() const {
            assert(!empty());

            return id_vector[0];
        }
    };

    class TrieNode {
    public:
        using pointer_type = TrieNode *;
        static const size_t INF = std::numeric_limits<size_t>::max();
        std::array<pointer_type, card> children;

        TrieNode() : target_node{nullptr}, ids{nullptr} {
            children.fill(nullptr);
        }

        pointer_type target_node;

        IdCounter *ids;

        void compress_to_shortest(pointer_type p = nullptr, size_t dist = 0) {
            if (!p && ids) {
                p = this;
                dist = 0;
            }

            if (ids) {
                target_node = p;
            }

            for (auto &child : children) {
                if (child) {
                    child->compress_to_shortest(p, dist + 1);
                }
            }
        }

        void compress_to_itselft() {

            if (ids) {
                target_node = this;
            }

            for (auto &child : children) {
                if (child) {
                    child->compress_to_itselft();
                }
            }
        }

        void checkout(std::unordered_map<size_t, size_t> &result) const {
            if (ids && !ids->empty()) {
                size_t id = target_node->ids->represent();
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
                assert(!target_node->ids->empty());
                size_t id = target_node->ids->represent();
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
};
}
// vim: ts=4:sw=4
