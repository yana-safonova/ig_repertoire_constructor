#pragma once

#include <vector>
#include <memory>
#include <utility>

#include <boost/unordered_map.hpp>
#include <seqan/sequence.h>

#include "pog_parameters.hpp"

namespace pog {

using pair_vector = std::vector<std::pair<size_t, size_t>>;

size_t hamming_distance(seq_t const& s1, seq_t const& s2);

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
    node(boost::unordered_map<size_t, std::shared_ptr<node>> const& id_map);
    node(kmer const& source, boost::unordered_map<size_t, std::shared_ptr<node>> const& id_map);

    node(node const&) = delete;
    node& operator=(node const&) = delete;

    void add_read();
    void set_coverage(float coverage);
    void add_output_edge(std::shared_ptr<node> const& next, float coverage = 1.f);
    void remove_output_edge(std::shared_ptr<node> const& next);
    bool on_upath();
    bool on_bulge();
    void remove_node();
    static bool join_nodes(std::shared_ptr<node>& a, std::shared_ptr<node>& b,
                           boost::unordered_map<size_t, std::shared_ptr<node>> const& id_map);

    bool dummy() const noexcept;
    bool equals(kmer const& potential_match) const noexcept;
    float coverage() const noexcept;
    size_t get_length() const noexcept;
    boost::unordered_map<size_t, float> const& get_input_edges() const noexcept;
    boost::unordered_map<size_t, float> const& get_output_edges() const noexcept;
    seq_t const& get_sequence() const noexcept;

    void set_index(size_t index);
    size_t get_index() const noexcept;
    size_t get_id() const noexcept;

private:

    static size_t global_id;

    seq_t sequence_;
    u64 hash_;
    float coverage_;
    size_t index_;
    size_t const id_;

    //                  node id, coverage
    boost::unordered_map<size_t, float> input_edges_;
    boost::unordered_map<size_t, float> output_edges_;
    boost::unordered_map<size_t, std::shared_ptr<node>> const& id_map_;

};

} // namespace pog
