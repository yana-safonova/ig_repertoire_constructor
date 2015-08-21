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
  private:
    static constexpr size_t card = seqan::ValueSize<Tletter>::VALUE;

    class TrieNode {
      public:
        using pointer_type = TrieNode*;
        static const size_t INFu = -1u;
        std::array<pointer_type, card> children;

        TrieNode() : nearest_leaf{nullptr}, nearest_leaf_distance{INFu}, ids{nullptr} {
          children.fill(nullptr);
        }

        pointer_type nearest_leaf;
        size_t nearest_leaf_distance;
        std::vector<size_t> *ids;

        void compute_nearest_leaf_distance() {
          nearest_leaf_distance = INFu;

          for (auto &child : children) {
            if (child) {
              child->compute_nearest_leaf_distance();
              if (nearest_leaf_distance > child->nearest_leaf_distance + 1) {
                nearest_leaf_distance = child->nearest_leaf_distance;
                nearest_leaf = child->nearest_leaf;
              }
            }
          }

          if (nearest_leaf_distance == INFu) {
            nearest_leaf_distance = 0;
            nearest_leaf = this;
          }
        }

        void checkout(std::unordered_map<size_t, size_t> &result) const {
          if (ids && !ids->empty()) {
            assert(!nearest_leaf->ids->empty());
            size_t id = nearest_leaf->ids->at(0);
            result[id] += ids->size();
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

  public:
    Trie() = default;
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;
    Trie(Trie&&) = default;
    Trie& operator=(Trie&&) = default;

    template<typename Tcont>
      Trie(const Tcont &cont) {
        size_t i = 0;
        for (const auto &s : cont) {
          add(s, i);
          ++i;
        }
      }

    template<typename T, typename Tf>
      void add(const T &s, size_t id, const Tf &toIndex) {
        if (!root) {
          root.reset(new TrieNode);
        }

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

    std::unordered_map<size_t, size_t> checkout(size_t nbucket) {
      root->compute_nearest_leaf_distance();

      std::unordered_map<size_t, size_t> result(nbucket);
      root->checkout(result);

      return result;
    }

    std::unordered_map<size_t, size_t> checkout() {
      size_t nleaves = root->leaves_count();
      return checkout(nleaves);
    }
};

// vim: ts=4:sw=4
