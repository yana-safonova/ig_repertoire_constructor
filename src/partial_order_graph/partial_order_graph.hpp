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
    
    void load(std::string const& prefix);
    void add_sequence(seq_t const& sequence, id_t const& read_id);
    seq_t correct_read(seq_t const& sequence, id_t const& read_id) const;
    void shrink_upaths();
    void shrink_bulges();
    void remove_low_covered();

    void save_dot(std::string const& filename, bool print_sequences = false) const;
    void save_nodes(std::string const& filename) const;
    size_t nodes_count() const noexcept;
    float total_coverage() const noexcept;

    void clear();

private:

    void align_cell(std::vector<kmer> const& kmer_seq, std::vector<float>& scores,
                    pair_vector& transitions, size_t i, size_t j, size_t m, bool relative_mismatch_penalty) const;
    size_t starting_point(std::vector<float> const& scores, size_t m) const;
    pair_vector select_matches(std::vector<float> const& scores, pair_vector const& transitions, size_t m) const;
    seq_t rebuild_read(size_t start, pair_vector const& transitions, size_t m) const;
    std::pair<std::vector<float>, pair_vector> align(std::vector<kmer> const& kmer_seq, bool relative_mismatch_penalty) const;

    void add_mismatching_region(std::vector<std::shared_ptr<node>>& new_nodes, std::vector<kmer> const& kmer_seq,
                                size_t i1, size_t i2, size_t j1, size_t j2);
    void update_nodes(std::vector<kmer> const& kmer_seq, pair_vector const& matches);
    void update_nodes_indexes();
    void clean_nodes();

    std::vector<size_t> most_covered_path() const;
    std::vector<size_t> restore_path(std::vector<size_t> const& next) const;

    std::vector<std::shared_ptr<node>> nodes_;
    boost::unordered_map<size_t, std::shared_ptr<node>> id_map_;

};

void from_file(partial_order_graph& graph, std::string const& filename, bool use_forward_reads, bool use_rc_reads);
void correct_reads(partial_order_graph const& graph, std::string const& filename_in, std::string const& filename_out,
                   bool use_forward_reads, bool use_rc_reads);
void save_graph(partial_order_graph const& graph, std::string const& prefix, bool show_sequences);
void draw_graph(std::string const& prefix);

} // namespace pog
