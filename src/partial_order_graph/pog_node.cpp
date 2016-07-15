#include "pog_node.hpp"

namespace pog {

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

void node::on_upath() {
    if (input_edges_.size() != 1 || dummy())
        return;

    auto const& input_edge = *input_edges_.begin();
    node* prev = input_edge.first;
    if (prev->output_edges_.size() != 1 || prev->dummy())
        return;

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
}

bool node::join_with_brother(node* brother) {
    if (coverage_ < pog_parameters::instance().bulge_coverage_difference * brother->coverage_ ||
            get_length() < brother->get_length())
        return false;
    TRACE("Shrinking bulge. Coverages: " << coverage_ << " and " << brother->coverage_
            << "; lengths: " << get_length() << " and " << brother->get_length());
    coverage_ += brother->coverage_;
    brother->sequence_ = "";
    for (auto& entry : brother->input_edges_)
        entry.first->output_edges_.erase(brother);
    for (auto& entry : brother->output_edges_)
        entry.first->input_edges_.erase(brother);
    return true;
}

void node::on_bulge() {
    if (dummy() || input_edges_.size() != 1 || output_edges_.size() != 1)
        return;

    auto const& input_edge = *input_edges_.begin();
    node* prev = input_edge.first;
    auto const& output_edge = *output_edges_.begin();
    node* next = output_edge.first;

    for (auto const& entry : prev->output_edges_) {
        node* brother = entry.first;
        if (brother == this || brother->output_edges_.find(next) == brother->output_edges_.end())
            continue;
        if (brother->input_edges_.size() != 1 || brother->output_edges_.size() != 1)
            continue;
        if (join_with_brother(brother)) {
            on_bulge();
            on_upath();
            next->on_upath();
            return;
        } else if (brother->join_with_brother(this)) {
            brother->on_bulge();
            brother->on_upath();
            next->on_upath();
            return;
        }
    }
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
