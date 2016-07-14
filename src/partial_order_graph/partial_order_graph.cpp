#include "partial_order_graph.hpp"

#include <fstream>

namespace pog {

partial_order_graph::partial_order_graph() {
    nodes_.push_back(new node());
    nodes_.push_back(new node());
    update_nodes_indexes();
}

partial_order_graph::~partial_order_graph() {
    for (node* v : nodes_)
        delete v;
}

void partial_order_graph::add_sequence(seq_t const& sequence, id_t const& read_id) {
    TRACE("Adding " << read_id << " to the graph");
    if (length(sequence) < pog_parameters::instance().get_kmer_size()) {
        DEBUG(read_id << " is too short: " << length(sequence));
        return;
    }

    std::vector<kmer> kmer_seq = sequence_to_kmers(sequence);
    TRACE("\tThere are " << kmer_seq.size() << " consecutive kmers in the sequence");

    pair_vector matches = align(kmer_seq);
    update_nodes(kmer_seq, matches);
    update_nodes_indexes();
}

void partial_order_graph::align_cell(std::vector<kmer> const& kmer_seq, std::vector<float>& scores,
                                     pair_vector& transitions, size_t i, size_t j, size_t m) const {
    static pog_parameters& parameters = pog_parameters::instance();

    float align_score = 0;
    if (j < m)
        align_score = nodes_[i + 1]->equals(kmer_seq[j]) ? 1 : parameters.mismatch_penalty;

    float& current_score = scores[i * (m + 1) + j];
    current_score = -1e6f;
    std::pair<size_t, size_t>& current_transition = transitions[i * (m + 1) + j];

    if (j < m && scores[i * (m + 1) + j + 1] + parameters.gap_penalty > current_score) {
        current_score = scores[i * (m + 1) + j + 1] + parameters.gap_penalty;
        current_transition = std::make_pair(i, j + 1);
    }
    for (auto const& output_edge : nodes_[i + 1]->get_output_edges()) {
        size_t k = node_indexes_.at(output_edge.first) - 1;
        if (j < m && scores[k * (m + 1) + j + 1] + align_score > current_score) {
            current_score = scores[k * (m + 1) + j + 1] + align_score;
            current_transition = std::make_pair(k, j + 1);
        }
        if (scores[k * (m + 1) + j] + parameters.gap_penalty > current_score) {
            current_score = scores[k * (m + 1) + j] + parameters.gap_penalty;
            current_transition = std::make_pair(k, j);
        }
    }
}

pair_vector partial_order_graph::select_matches(std::vector<float> const& scores,
                                                pair_vector const& transitions, size_t m) const {
    size_t n = nodes_.size() - 2;
    size_t i = 0;
    size_t j = 0;
    pair_vector matches;

    while (i != n || j != m) {
        size_t i_next;
        size_t j_next;
        std::tie(i_next, j_next) = transitions[i * (m + 1) + j];
        if (scores[i * (m + 1) + j] == scores[i_next * (m + 1) + j_next] + 1) {
            matches.push_back(std::make_pair(i, j));
        }
        i = i_next;
        j = j_next;
    }
    TRACE("\tMatches: " << matches.size());
    return matches;
}

// Return value: Vector of matches + mismatches (i, j), i in graph, j in kmer_seq
pair_vector partial_order_graph::align(std::vector<kmer> const& kmer_seq) const {
    size_t n = nodes_.size() - 2;
    size_t m = kmer_seq.size();

    std::vector<float> scores((n + 1) * (m + 1));
    pair_vector transitions((n + 1) * (m + 1));

    scores[n * (m + 1) + m] = 0;
    for (size_t i = n; i <= n; --i) {
        for (size_t j = m; j <= m; --j) {
            if (i < n || j < m)
                align_cell(kmer_seq, scores, transitions, i, j, m);
        }
    }

    TRACE("\tScore: " << scores[0]);
    return select_matches(scores, transitions, m);
}

// Inserting nodes [i1, i2) from graph, then [j1, j2) from kmer_seq, then i2 from graph.
// Increasing coverage of every node, except i2
// If j1 >= j2 adding edge between i1, i2
void partial_order_graph::add_mismatching_region(std::vector<node*>& new_nodes,
            std::vector<kmer> const& kmer_seq, size_t i1, size_t i2, size_t j1, size_t j2) {
    for (size_t i = i1; i < i2; ++i)
        new_nodes.push_back(nodes_[i]);
    nodes_[i1]->add_read();

    if (j1 < j2) {
        new_nodes.push_back(new node(kmer_seq[j1]));
        nodes_[i1]->add_output_edge(new_nodes.back());
    }

    for (size_t j = j1 + 1; j < j2; ++j) {
        node* next = new node(kmer_seq[j]);
        new_nodes.back()->add_output_edge(next);
        new_nodes.push_back(next);
    }
    new_nodes.back()->add_output_edge(nodes_[i2]);
}

void partial_order_graph::update_nodes(std::vector<kmer> const& kmer_seq,
                                       pair_vector const& matches) {
    size_t n = nodes_.size() - 1;
    size_t m = kmer_seq.size();
    std::vector<node*> new_nodes;

    size_t i_prev = 0;
    size_t j_prev = 0;
    for (auto const& match : matches) {
        size_t i = match.first + 1;
        size_t j = match.second;
        add_mismatching_region(new_nodes, kmer_seq, i_prev, i, j_prev, j);
        i_prev = i;
        j_prev = j + 1;
    }

    add_mismatching_region(new_nodes, kmer_seq, i_prev, n, j_prev, m);
    new_nodes.push_back(nodes_.back());
    new_nodes.back()->add_read();
    nodes_ = new_nodes;
}

void partial_order_graph::update_nodes_indexes() {
    DEBUG("Nodes: " << nodes_.size());
    node_indexes_.clear();
    for (size_t i = 0; i < nodes_.size(); ++i) {
        node_indexes_[nodes_[i]] = i;
    }
}

void partial_order_graph::save_dot(std::string const& filename, std::string const& graph_name, bool print_sequences) const {
    std::ofstream dot_file(filename);
    if (!dot_file) {
        ERROR("Could not write to " << filename);
        return;
    }
    dot_file << "digraph " << graph_name << " {\n";
    for (size_t i = 0; i < nodes_.size(); ++i) {
        dot_file << '\t' << i << " [label = \"#" << i << "\\n" << nodes_[i]->coverage();
        if (print_sequences && length(nodes_[i]->get_sequence()))
            dot_file << ": " << nodes_[i]->get_sequence();
        dot_file << "\"]\n";

        for (auto const& entry : nodes_[i]->get_output_edges()) {
            dot_file << '\t' << i << " -> " << node_indexes_.at(entry.first) << " [label = \"" << entry.second << "\"]\n";
        }
    }
    dot_file <<"}\n";
}

void partial_order_graph::save_nodes(std::string const& filename) const {
    std::ofstream nodes_file(filename);
    if (!nodes_file) {
        ERROR("Could not write to " << filename);
        return;
    }

    nodes_file << "node\tsequence\tcoverage\n";
    for (size_t i = 1; i < nodes_.size() - 1; ++i)
        nodes_file << i << '\t' << nodes_[i]->get_sequence() << '\t' << nodes_[i]->coverage() << '\n';
}

size_t partial_order_graph::nodes_count() const noexcept {
    return nodes_.size();
}

partial_order_graph from_file(std::string const& filename, bool use_reverse_complement) {
    seqan::SeqFileIn reader(filename.c_str());

    if (!guessFormat(reader)) {
        ERROR("Could not detect file format " << filename);
        exit(1);
    }

    id_t id;
    seq_t seq;
    if (atEnd(reader)) {
        ERROR("Cannot read from file " << filename);
        exit(1);
    }

    INFO("Starting reading from " << filename);

    partial_order_graph graph;
    size_t count = 0;
    while (!atEnd(reader))
    {
        readRecord(id, seq, reader);
        DEBUG(++count << ": " << id);
        graph.add_sequence(seq, id);

        if (use_reverse_complement) {
            reverseComplement(seq);
            graph.add_sequence(seq, id);
        }
    }

    INFO(filename << ": " << count << " reads");
    return graph;
}

} // namespace pog
