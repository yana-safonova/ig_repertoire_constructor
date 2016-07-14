#pragma once

#include <vector>
#include <utility>
#include <functional>

#include <boost/unordered_map.hpp>
#include <seqan/sequence.h>

#include "pog_parameters.hpp"

namespace pog {

using pair_vector = std::vector<std::pair<size_t, size_t>>;

struct kmer {
    // First k letters
    kmer(seq_t const& read_sequence);
    kmer(seq_t const& read_sequence, size_t start, u64 hash) noexcept;
    nt_t operator[](size_t i) const;
    u64 get_hash() const noexcept;
    size_t get_start() const noexcept;

    seq_t const& read_sequence;

private:

    size_t start_;
    u64 hash_;

};

std::vector<kmer> sequence_to_kmers(seq_t const& sequence);

struct subnode {
    subnode(kmer const& source, size_t read_number);

    void add_read(size_t read_number, size_t position);
    pair_vector const& get_reads() const noexcept;
    size_t coverage() const noexcept;
    seq_t const& get_sequence() const noexcept;
    bool equals(kmer const& other) const;

private:

    seq_t sequence_;
    u64 hash_;
    pair_vector reads_;
};

struct node {

    node();
    node(kmer const& source, size_t read_number);
    node(node const&) = delete;
    node& operator=(node const&) = delete;

    void add_kmer(kmer const& source, size_t read_number);
    void add_output_edge(node* next);
    bool contains(kmer const& potential_match) const;
    boost::unordered_map<node*, size_t> const& get_output_edges() const noexcept;
    boost::unordered_map<u64, std::vector<subnode>> const& get_subnodes() const noexcept;
    void for_every_subnode(std::function<void(subnode const&)> f) const;

private:

    boost::unordered_map<u64, std::vector<subnode>> subnodes_;
    boost::unordered_map<node*, size_t> input_edges_;
    boost::unordered_map<node*, size_t> output_edges_;

};

} // namespace pog
