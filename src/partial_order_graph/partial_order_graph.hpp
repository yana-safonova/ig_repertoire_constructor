#pragma once

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <string>
#include <boost/unordered_map.hpp>
#include "pog_node.hpp"

namespace pog {


struct partial_order_graph {

    partial_order_graph();
    ~partial_order_graph();
    void add_sequence(seq_t const& sequence, id_t const& read_id);
    void shrink_upaths();
    void shrink_bulges();

    void save_dot(std::string const& filename, std::string const& graph_name,
                  bool print_sequences = false) const;
    void save_nodes(std::string const& filename) const;
    size_t nodes_count() const noexcept;

private:

    void align_cell(std::vector<kmer> const& kmer_seq, std::vector<float>& scores,
                    pair_vector& transitions, size_t i, size_t j, size_t m) const;
    pair_vector select_matches(std::vector<float> const& scores, pair_vector const& transitions, size_t m) const;
    pair_vector align(std::vector<kmer> const& kmer_seq) const;

    void add_mismatching_region(std::vector<node*>& new_nodes, std::vector<kmer> const& kmer_seq,
                                size_t i1, size_t i2, size_t j1, size_t j2);
    void update_nodes(std::vector<kmer> const& kmer_seq, pair_vector const& matches);
    void update_nodes_indexes();
    void clean_nodes();
    std::vector<size_t> most_covered_path() const;

    std::vector<node*> nodes_;
    boost::unordered_map<node*, size_t> node_indexes_;

};

partial_order_graph from_file(std::string const& filename, bool use_reverse_complement = true);

} // namespace pog
