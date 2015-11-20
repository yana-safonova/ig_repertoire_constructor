#pragma once

#include <vector>
#include <array>
#include <memory>
#include <cassert>
#include <algorithm>
#include <unordered_map>

#include <seqan/seq_io.h>
using seqan::length;

template<typename Tletter = seqan::Dna5>
class Trie {
public:
    Trie() : root{new TrieNode} {  }
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;
    Trie(Trie&&) = default;
    Trie& operator=(Trie&&) = default;
    ~Trie() = default;

    template<typename Tcont>
    Trie(const Tcont &cont) {
        root.reset(new TrieNode);
        size_t i = 0;
        for (const auto &s : cont) {
            add(s, i);
            ++i;
        }
    }

    template<typename T, typename Tf>
    void add(const T &s, size_t id, const Tf &toIndex) {
        assert(!isCompressed());

        typename TrieNode::pointer_type p = this->root.get();

        for (size_t i = 0; i < length(s); ++i) {
            size_t el = toIndex(s[i]);
            assert((0 <= el) && (el < p->children.size()));

            if (!p->children[el]) {
                p->children[el] = new TrieNode();
            }

            p = p->children[el];
        }

        if (!p->ids) {
            p->ids = new std::vector<size_t>;
        }

        p->ids->push_back(id);
    }

    template<typename T>
    void add(const T &s, size_t id) {
        auto to_size_t = [](const Tletter &letter) -> size_t { return seqan::ordValue(letter); };
        add(s, id, to_size_t);
    }

    bool isCompressed() const {
        return compressed;
    }

    void compress() {
        if (!isCompressed()) {
            // root->compress_to_longest();
            root->compress_to_shortest();
        }
    }

    std::unordered_map<size_t, size_t> checkout(size_t nbucket = 0) {
        if (nbucket == 0) {
            nbucket = root->leaves_count();
        }

        if (!isCompressed()) compress();

        std::unordered_map<size_t, size_t> result(nbucket);
        root->checkout(result);

        return result;
    }

    std::unordered_map<size_t, std::vector<size_t>> checkout_ids(size_t nbucket = 0) {
        if (nbucket == 0) {
            nbucket = root->leaves_count();
        }

        if (!isCompressed()) compress();

        std::unordered_map<size_t, std::vector<size_t>> result(nbucket);
        root->checkout(result);

        return result;
    }

private:
    static constexpr size_t card = seqan::ValueSize<Tletter>::VALUE;

    class TrieNode {
    public:
        using pointer_type = TrieNode*;
        static const size_t INFu = -1u;
        std::array<pointer_type, card> children;

        TrieNode() : target_node{nullptr}, target_node_distance{INFu}, ids{nullptr} {
            children.fill(nullptr);
        }

        pointer_type target_node;

        size_t target_node_distance;
        std::vector<size_t> *ids;

        void compress_to_longest() {
            target_node_distance = INFu;

            for (auto &child : children) {
                if (child) {
                    child->compress_to_longest();
                    if (target_node_distance > child->target_node_distance + 1) {
                        target_node_distance = child->target_node_distance;
                        target_node = child->target_node;
                    }
                }
            }

            if (target_node_distance == INFu) {
                target_node_distance = 0;
                target_node = this;
            }
        }

        void compress_to_shortest(pointer_type p = nullptr, size_t dist = 0) {
            target_node_distance = INFu;

            if (!p && ids) {
                p = this;
                dist = 0;
            }

            if (ids) {
                target_node = p;
                target_node_distance = dist;
            }

            for (auto &child : children) {
                if (child) {
                    child->compress_to_shortest(p, dist + 1);
                }
            }
        }

        void compress_to_itselft() {
            target_node_distance = INFu;

            if (ids) {
                target_node = this;
                target_node_distance = 0;
            }

            for (auto &child : children) {
                if (child) {
                    child->compress_to_itselft();
                }
            }
        }

        void checkout(std::unordered_map<size_t, size_t> &result) const {
            if (ids && !ids->empty()) {
                assert(!target_node->ids->empty());
                size_t id = target_node->ids->at(0);
                result[id] += ids->size();
            }

            // DFS
            for (const auto &child : children) {
                if (child) child->checkout(result);
            }
        }

        void checkout(std::unordered_map<size_t, std::vector<size_t>> &result) const {
            if (ids && !ids->empty()) {
                assert(!target_node->ids->empty());
                size_t id = target_node->ids->at(0);
                result[id].insert(result[id].end(), ids->cbegin(), ids->cend());
            }

            // DFS
            for (const auto &child : children) {
                if (child) child->checkout(result);
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

            for (auto &child : children) {
                if (child) delete child;
            }
        }
    };

    std::unique_ptr<TrieNode> root;
    bool compressed = false;
};

// vim: ts=4:sw=4
