#include "pog_node.hpp"

namespace pog {

size_t hamming_distance(seq_t const& s1, seq_t const& s2) {
    if (length(s1) != length(s2)) {
        WARN("Different lengths, cannot compute Hamming distance: " << length(s1) << " " << length(s2));
        return -1;
    }
    size_t distance = 0;
    for (size_t i = 0; i < length(s1); ++i)
        distance += s1[i] != s2[i] ? 1 : 0;
    return distance;
}

// ----------------- kmer

kmer::kmer(seq_t const& read_sequence)
        : read_sequence(read_sequence)
        , start_(0)
        , hash_(0) {
    assert(pog_parameters::instance().get_kmer_size() <= seqan::length(read_sequence));
    for (size_t i = 0; i < pog_parameters::instance().get_kmer_size(); ++i) {
        hash_ = hash_ * seqan::ValueSize<nt_t>::VALUE // i.e. count of nt types (A,C,G,T,N)
                + ordValue(read_sequence[i]);
    }
}

kmer::kmer(seq_t const& read_sequence, size_t start, u64 hash) noexcept
        : read_sequence(read_sequence)
        , start_(start)
        , hash_(hash) {
    // EMPTY
}

nt_t kmer::operator[](size_t i) const {
    assert(start_ + i < seqan::length(read_sequence));
    return read_sequence[start_ + i];
}

u64 kmer::get_hash() const noexcept {
    return hash_;
}

size_t kmer::get_start() const noexcept {
    return start_;
}

std::vector<kmer> sequence_to_kmers(seq_t const& sequence) {
    std::vector<kmer> kmers;
    pog_parameters& parameters = pog_parameters::instance();
    if (seqan::length(sequence) < parameters.get_kmer_size())
        return kmers;

    kmers.push_back(kmer(sequence));
    u64 hash = kmers.back().get_hash();
    for (size_t i = parameters.get_kmer_size(); i < seqan::length(sequence); ++i) {
        hash = (hash * parameters.alphabet_size + ordValue(sequence[i])) % parameters.mask;
        kmers.push_back(kmer(sequence, i - parameters.get_kmer_size() + 1, hash));
    }
    return kmers;
}


// ----------------- node

node::node()
        : hash_(-1)
        , coverage_(0) {
    // EMPTY
}

node::node(kmer const& source)
        : sequence_(seqan::infixWithLength(source.read_sequence,
                                           source.get_start(),
                                           pog_parameters::instance().get_kmer_size()))
        , hash_(source.get_hash())
        , coverage_(1) {
    // EMPTY
}

void node::add_read() {
    ++coverage_;
}

void node::add_output_edge(node* next, float coverage) {
    output_edges_[next] += coverage;
    next->input_edges_[this] += coverage;
}

bool node::on_upath() {
    if (input_edges_.size() != 1 || dummy())
        return false;

    auto const& input_edge = *input_edges_.begin();
    node* prev = input_edge.first;
    if (prev->output_edges_.size() != 1 || prev->dummy())
        return false;

    if (prev->coverage_ != coverage_) {
        DEBUG("Coverage is different: " << prev->coverage_ << " != " << coverage_);
    }

    float prev_length = static_cast<float>(prev->get_length());
    float curr_length = static_cast<float>(get_length());
    prev->coverage_ = (prev->coverage_ * prev_length + coverage_ * curr_length) / (prev_length + curr_length);

    append(prev->sequence_, suffix(sequence_, pog_parameters::instance().get_kmer_size() - 1));
    sequence_ = "";

    for (auto const& entry : output_edges_) {
        prev->add_output_edge(entry.first, std::min(input_edge.second, entry.second));
        entry.first->input_edges_.erase(this);
    }
    prev->output_edges_.erase(this);
    return true;
}

bool node::join_nodes(node* a, node* b) {
    if (a->get_length() != b->get_length() || a->dummy())
        return false;
    if (a->output_edges_.find(b) != a->output_edges_.end() || b->output_edges_.find(a) != b->output_edges_.end())
        return false;
    TRACE("\tHamming: " << hamming_distance(a->sequence_, b->sequence_));
    TRACE("\tMargin:  " << std::ceil(static_cast<float>(length(a->sequence_)) * pog_parameters::instance().bulges_hamming_ratio));
    if (static_cast<float>(hamming_distance(a->sequence_, b->sequence_))
            > std::ceil(static_cast<float>(length(a->sequence_)) * pog_parameters::instance().bulges_hamming_ratio))
        return false;

    if (a->input_edges_.size() + a->output_edges_.size() < b->input_edges_.size() + b->output_edges_.size())
        std::swap(a, b);
    TRACE("Joining nodes. Coverage: " << a->coverage_ << " and " << b->coverage_);

    for (auto& input_edge : b->input_edges_) {
        input_edge.first->output_edges_.erase(b);
        input_edge.first->add_output_edge(a, input_edge.second);
    }
    for (auto& output_edge : b->output_edges_) {
        output_edge.first->input_edges_.erase(b);
        a->add_output_edge(output_edge.first, output_edge.second);
    }
    a->coverage_ += b->coverage_;
    b->sequence_ = "";
    return true;
}

bool process_brothers(std::vector<node*>& brothers) {
    size_t n = brothers.size();
    bool flag = false;
    for (size_t i = 0; i < n - 1; ++i)
        for (size_t j = i + 1; j < n; ++j)
            flag = node::join_nodes(brothers[i], brothers[j]) || flag;
    return flag;
}

void node::on_bulge() {
    //            node via one, intermediate nodes
    boost::unordered_map<node*, std::vector<node*>> intermediate_nodes;
    for (auto const& entry1 : output_edges_)
        for (auto const& entry2 : entry1.first->output_edges_)
            intermediate_nodes[entry2.first].push_back(entry1.first);

    for (auto& entry : intermediate_nodes) {
        if (entry.second.size() > 1 && process_brothers(entry.second)) {
            entry.first->on_upath();
        }
    }

    for (auto const& entry : output_edges_)
        entry.first->on_upath();
}

bool node::dummy() const noexcept {
    return !length(sequence_);
}

bool node::equals(kmer const& potential_match) const noexcept {
    if (hash_ != potential_match.get_hash() || dummy())
        return false;

    assert(length(sequence_) == pog_parameters::instance().get_kmer_size());
    for (size_t i = 0; i < pog_parameters::instance().get_kmer_size(); ++i) {
        if (sequence_[i] != potential_match[i])
            return false;
    }
    return true;
}

float node::coverage() const noexcept {
    return coverage_;
}

size_t node::get_length() const noexcept {
    return length(sequence_) - pog_parameters::instance().get_kmer_size() + 1;
}

boost::unordered_map<node*, float> const& node::get_input_edges() const noexcept {
    return input_edges_;
}

boost::unordered_map<node*, float> const& node::get_output_edges() const noexcept {
    return output_edges_;
}

seq_t const& node::get_sequence() const noexcept {
    return sequence_;
}

} // namespace pog
