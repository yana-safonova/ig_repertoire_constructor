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

node::node() {
    hash_ = static_cast<u64>(-1);
}

node::node(kmer const& source, directed_seqan_read* read)
        : sequence_(seqan::infixWithLength(source.read_sequence,
                                           source.get_start(),
                                           pog_parameters::instance().get_kmer_size())) {
    add_read(read, source.get_start());
    hash_ = source.get_hash();
}

void node::add_read(directed_seqan_read* read, size_t position) {
    reads_.push_back(std::make_pair(read, position));
}

void node::add_output_edge(node* next) {
    auto it = output_edges_.find(next);
    if (it == output_edges_.end()) {
        output_edges_[next] = 1;
        next->input_edges_[this] = 1;
    } else {
        output_edges_[next] += 1;
        next->input_edges_[this] += 1;
    }
}

bool node::dummy() const noexcept {
    return !length(sequence_);
}

bool node::sequences_equal(kmer const& potential_match) const noexcept {
    if (hash_ != potential_match.get_hash() || dummy())
        return false;

    assert(length(sequence_) == pog_parameters::instance().get_kmer_size());
    for (size_t i = 0; i < pog_parameters::instance().get_kmer_size(); ++i) {
        if (sequence_[i] != potential_match[i])
            return false;
    }
    return true;
}

size_t node::coverage() const noexcept {
    return reads_.size();
}

std::vector<std::pair<directed_seqan_read*, size_t>> const& node::get_reads() const noexcept {
    return reads_;
}

boost::unordered_map<node*, size_t> const& node::get_output_edges() const noexcept {
    return output_edges_;
}

seq_t const& node::get_sequence() const noexcept {
    return sequence_;
}

} // namespace pog
