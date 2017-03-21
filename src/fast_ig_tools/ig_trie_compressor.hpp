#pragma once

#include <algorithm>
#include <array>
#include <boost/pool/object_pool.hpp>
#include <cassert>
#include <limits>
#include <memory>
#include <boost/unordered_map.hpp>
#include <vector>

#include <seqan/seq_io.h>

namespace fast_ig_tools {
template <typename T>
using Decay = typename std::decay<T>::type;

template <typename TArray>
using ValueType = Decay<decltype(Decay<TArray>()[0])>;

class Compressor {
public:
    enum class Type {HashCompressor, TrieCompressor};
    virtual std::vector<size_t> checkout() = 0;
    virtual ~Compressor() = default;

    template <typename TValue, typename... Args>
    static std::unique_ptr<Compressor> factor(Compressor::Type type, Args&&... args);

    template <typename TIter>
    static std::vector<size_t> compressed_reads_indices(TIter b, TIter e, Compressor::Type type = Compressor::Type::TrieCompressor) {
		using TValue = ValueType<decltype(*(TIter()))>;
        auto compressor = Compressor::factor<TValue>(type, b, e);

        return compressor->checkout();
    }

    template <typename TVector>
    static std::vector<size_t> compressed_reads_indices(const TVector &reads, Compressor::Type type = Compressor::Type::TrieCompressor) {
        return compressed_reads_indices(reads.cbegin(), reads.cend(), type);
    }

    template <typename TIter>
    static auto compressed_reads(TIter b, TIter e, Compressor::Type type = Compressor::Type::TrieCompressor) -> std::vector<Decay<decltype(*(TIter()))>>  {
        auto indices = compressed_reads_indices(b, e, type);

        std::vector<Decay<decltype(*(TIter()))>> result;
        for (size_t i = 0; i < indices.size(); ++i, ++b) {
            if (indices[i] == i) {
                result.push_back(*b);
            }
        }

        return result;
    }

    template <typename TVector>
    static auto compressed_reads(const TVector &reads, Compressor::Type type = Compressor::Type::TrieCompressor) -> std::vector<ValueType<TVector>>  {
        return compressed_reads(reads.cbegin(), reads.cend(), type);
    }
};


class HashCompressor : public Compressor {
public:
    HashCompressor() {}
    HashCompressor(const HashCompressor &) = delete;
    HashCompressor &operator=(const HashCompressor &) = delete;
    HashCompressor(HashCompressor &&) = default;
    HashCompressor &operator=(HashCompressor &&) = default;
    virtual ~HashCompressor() = default;

    size_t size() const {
        return size_;
    }

    template <typename TCont>
    HashCompressor(const TCont &cont) : HashCompressor(cont.cbegin(), cont.cend()) { }

    template <typename TIter>
    HashCompressor(TIter b, TIter e) : HashCompressor() {
        for (; b !=e; ++b) {
            add(*b);
        }
    }
    template <typename T>
    void add(const T &s) {
        seqan::String<char> str = s;
        map_[seqan::toCString(str)].push_back(size_++);
    }

    virtual std::vector<size_t> checkout() {
        std::vector<size_t> result(this->size());
        for (const auto &kv : map_) {
            const auto &ids = kv.second;
            size_t target = ids[0];
            for (size_t id : ids) {
                result[id] = target;
            }
        }

        return result;
    }

private:
    boost::unordered_map<std::string, std::vector<size_t>> map_;
    size_t size_ = 0;
};


template <typename TValue = seqan::Dna5>
class TrieCompressor : public Compressor {
public:
    TrieCompressor() {
        // Do not use :member initialization syntax here. pool should be initialized BEFORE root
        root_ = pool_.construct();
    }
    TrieCompressor(const TrieCompressor &) = delete;
    TrieCompressor &operator=(const TrieCompressor &) = delete;
    TrieCompressor(TrieCompressor &&) = default;
    TrieCompressor &operator=(TrieCompressor &&) = default;
    virtual ~TrieCompressor() = default;


    size_t size() const {
        return size_;
    }

    template <typename TCont>
    TrieCompressor(const TCont &cont) : TrieCompressor(cont.cbegin(), cont.cend()) { }

    template <typename TIter>
    TrieCompressor(TIter b, TIter e) : TrieCompressor() {
        for (; b !=e; ++b) {
            add(*b);
        }
    }

    template <typename T>
    void add(const T &s) {
        assert(!isCompressed());

        TrieNode *p = root_;

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
                p->children[el] = pool_.construct();
            }

            p = p->children[el];
        }

        if (!p->ids) {
            p->ids = new IdCounter;
        }

        p->ids->add(size_);
        size_++;
    }

    bool isCompressed() const {
        return compressed_;
    }

    void compress() {
        if (!isCompressed()) {
            root_->compress_to_prefix();
        }
    }

    virtual std::vector<size_t> checkout() {
        if (!isCompressed()) {
            compress();
        }

        std::vector<size_t> result(this->size());
        root_->checkout_vector(result);

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

    boost::object_pool<TrieNode> pool_;
    TrieNode *root_;
    bool compressed_ = false;
    size_t size_ = 0;
};


template <typename TValue, typename... Args>
std::unique_ptr<Compressor> Compressor::factor(Compressor::Type type, Args&&... args) {
    switch (type) {
        case Compressor::Type::HashCompressor:
            return std::unique_ptr<Compressor>(new HashCompressor(std::forward<Args>(args)...));
        case Compressor::Type::TrieCompressor:
            return std::unique_ptr<Compressor>(new TrieCompressor<TValue>(std::forward<Args>(args)...));
        default:
            return std::unique_ptr<Compressor>(nullptr);
    }
}
}
// vim: ts=4:sw=4
