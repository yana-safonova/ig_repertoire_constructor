#pragma once

#include <vector>
#include <utility>

#include <boost/unordered_map.hpp>
#include <seqan/sequence.h>

#include "seqan_read.hpp"
#include "pog_parameters.hpp"

namespace pog {

struct node;

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

struct node {
    node();
    node(kmer const& source, directed_seqan_read* read);

    node(node const&) = delete;
    node& operator=(node const&) = delete;

    void add_read(directed_seqan_read* read, size_t position);
    void add_output_edge(node* next);

    bool dummy() const noexcept;
    bool sequences_equal(kmer const& potential_match) const noexcept;
    size_t coverage() const noexcept;
    std::vector<std::pair<directed_seqan_read*, size_t>> const& get_reads() const noexcept;
    boost::unordered_map<node*, size_t> const& get_output_edges() const noexcept;
    seq_t const& get_sequence() const noexcept;

private:

    seq_t sequence_;
    u64 hash_;

    //                                    read, sequence position in read
    std::vector<std::pair<directed_seqan_read*, size_t>> reads_;
    //                     node, coverage
    boost::unordered_map<node*, size_t> input_edges_;
    boost::unordered_map<node*, size_t> output_edges_;
};

} // namespace pog
