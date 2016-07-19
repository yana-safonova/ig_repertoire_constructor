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

void partial_order_graph::shrink_upaths() {
    for (size_t i = 1; i < nodes_.size() - 1; ++i) {
        nodes_[i]->on_upath();
    }
    clean_nodes();
    INFO(nodes_.size() << " nodes after shrinking unambigous paths");
}

void partial_order_graph::shrink_bulges() {
    bool flag = true;
    while (flag) {
        flag = false;
        for (size_t i = 1; i < nodes_.size(); ++i) {
            flag = nodes_[i]->on_bulge() || flag;
        }
        clean_nodes();
    }
    INFO(nodes_.size() << " nodes after shrinking bulges");
}

void partial_order_graph::remove_low_covered() {
    auto const& parameters = pog_parameters::instance();
    for (size_t i = 1; i < nodes_.size() - 1; ++i) {
        if (nodes_[i]->coverage() <= parameters.coverage_threshold)
            nodes_[i]->remove_node();

        std::vector<node*> to_remove;
        for (auto const& entry : nodes_[i]->get_output_edges()) {
            if (entry.second <= parameters.coverage_threshold)
                to_remove.push_back(entry.first);
        }
        for (node* v : to_remove)
            nodes_[i]->remove_output_edge(v);
    }
    clean_nodes();
    shrink_upaths();
    shrink_bulges();
    INFO(nodes_.size() << " nodes after removing low covered nodes and edges");
}

void partial_order_graph::clean_nodes() {
    std::vector<node*> new_nodes;
    new_nodes.push_back(nodes_[0]);

    for (size_t i = 1; i < nodes_.size() - 1; ++i) {
        if (nodes_[i]->dummy())
            delete nodes_[i];
        else
            new_nodes.push_back(nodes_[i]);
    }

    new_nodes.push_back(nodes_.back());
    nodes_ = new_nodes;
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
    float max_score = -1e6;
    size_t i = 0;

    for (auto const& entry : nodes_[0]->get_output_edges()) {
        size_t k = node_indexes_.at(entry.first) - 1;
        if (scores[k * (m + 1)] > max_score) {
            max_score = scores[k * (m + 1)];
            i = k;
        }
    }
    TRACE("\tScore: " << scores[i * (m + 1)]);
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

// Return value: Vector of matches (i, j), i in graph, j in kmer_seq
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
    if (j1 < j2)
        new_nodes.back()->add_output_edge(nodes_[i2]);
    else
        nodes_[i1]->add_output_edge(nodes_[i2]);
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

std::vector<size_t> partial_order_graph::most_covered_path() const {
    size_t n = nodes_.size();

    std::vector<float> coverage(n, 0);
    std::vector<size_t> next(n);

    for (size_t i = n - 1; i > 0; --i) {
        for (auto const& entry : nodes_[i]->get_input_edges()) {
            node* prev = entry.first;
            size_t k = node_indexes_.at(prev);
            if (coverage[k] < coverage[i] + entry.second) {
                coverage[k] = coverage[i] + entry.second;
                next[k] = i;
            }
        }
    }

    return restore_path(next);
}


std::vector<size_t> partial_order_graph::restore_path(std::vector<size_t> const& next) const {
    size_t n = nodes_.size();

    size_t sum_length = 0;
    float min_coverage = 1e10;
    float sum_coverage = 0;

    std::vector<size_t> path;
    path.push_back(0);
    path.push_back(next[0]);
    while (path.back() != n - 1) {
        min_coverage = std::min(min_coverage, nodes_[path.back()]->coverage());
        sum_length += nodes_[path.back()]->get_length();
        sum_coverage += nodes_[path.back()]->coverage();
        path.push_back(next[path.back()]);
    }
    INFO("Most covered path.");
    INFO("\tAverage coverage: " << sum_coverage / static_cast<float>(sum_length));
    INFO("\tSum length: " << sum_length);
    INFO("\tMin coverage on the path: " << min_coverage);

    return path;
}

void partial_order_graph::save_dot(std::string const& filename, bool print_sequences) const {
    std::ofstream dot_file(filename + ".dot");
    if (!dot_file) {
        ERROR("Could not write to " << filename << ".dot");
        return;
    }
    std::ofstream path_file(filename + ".fa");
    if (!path_file) {
        ERROR("Could not write to " << filename << ".fa");
        return;
    }
    path_file << ">Most covered path\n";

    std::vector<size_t> path = most_covered_path();
    auto it = path.begin();

    dot_file << "digraph {\n";
    for (size_t i = 0; i < nodes_.size(); ++i) {
        dot_file << '\t' << i << " [label=\"#" << i << "\\nCov: " << nodes_[i]->coverage();
        if (!nodes_[i]->dummy())
            dot_file << "\\nLen: " << nodes_[i]->get_length();
        if (print_sequences)
            dot_file << "\\n " << nodes_[i]->get_sequence();

        dot_file << "\"";
        if (*it == i) {
            dot_file << ", fillcolor=khaki1,style=filled";
            path_file << nodes_[i]->get_sequence();
            ++it;
        }
        dot_file << "]\n";

        for (auto const& entry : nodes_[i]->get_output_edges()) {
            dot_file << '\t' << i << " -> " << node_indexes_.at(entry.first) << " [label = \"" << entry.second << "\"]\n";
        }
    }
    dot_file <<"}\n";
    path_file << std::endl;
}

void partial_order_graph::save_nodes(std::string const& filename) const {
    std::ofstream nodes_file(filename);
    if (!nodes_file) {
        ERROR("Could not write to " << filename);
        return;
    }

    nodes_file << "node\tsequence\tlength\tcoverage\n";
    for (size_t i = 1; i < nodes_.size() - 1; ++i) {
        nodes_file << i << '\t' << nodes_[i]->get_sequence() << '\t' << nodes_[i]->get_length()
                   << '\t' << nodes_[i]->coverage() << '\n';
   }
}

size_t partial_order_graph::nodes_count() const noexcept {
    return nodes_.size();
}

partial_order_graph from_file(std::string const& filename, bool use_forward_reads, bool use_rc_reads) {
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
        ++count;
        DEBUG(count << ": " << id);
        if (use_forward_reads)
            graph.add_sequence(seq, id);

        if (use_rc_reads) {
            reverseComplement(seq);
            graph.add_sequence(seq, id);
        }
    }

    INFO(filename << ": " << count << " reads");
    return graph;
}

} // namespace pog
