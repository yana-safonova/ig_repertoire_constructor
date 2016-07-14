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

// ----------------- subnode

subnode::subnode(kmer const& source, size_t read_number)
        : sequence_(seqan::infixWithLength(source.read_sequence,
                                           source.get_start(),
                                           pog_parameters::instance().get_kmer_size()))
        , hash_(source.get_hash()) {
    add_read(read_number, source.get_start());
}

void subnode::add_read(size_t read_number, size_t position) {
    reads_.push_back(std::make_pair(read_number, position));
}

pair_vector const& subnode::get_reads() const noexcept {
    return reads_;
}

size_t subnode::coverage() const noexcept {
    return reads_.size();
}

seq_t const& subnode::get_sequence() const noexcept {
    return sequence_;
}

bool subnode::equals(kmer const& other) const {
    for (size_t i = 0; i < pog_parameters::instance().get_kmer_size(); ++i)
        if (sequence_[i] != other[i])
            return false;
    return true;
}

// ----------------- node

node::node() {
}

node::node(kmer const& source, size_t read_number) {
    add_kmer(source, read_number);
}

void node::add_kmer(kmer const& source, size_t read_number) {
    auto subnodes_it = subnodes_.find(source.get_hash());

    if (subnodes_it == subnodes_.end()) {
        subnodes_[source.get_hash()].push_back(subnode(source, read_number));
    } else {
        std::vector<subnode>& subnodes_vector = subnodes_it->second;
        for (subnode& v : subnodes_vector) {
            if (v.equals(source)) {
                v.add_read(read_number, source.get_start());
                return;
            }
        }
        subnodes_vector.push_back(subnode(source, read_number));

        DEBUG("Hashes are equal: hash = " << source.get_hash() << "; "
            << toCString(subnodes_vector.front().get_sequence()) << ", " << toCString(subnodes_vector.back().get_sequence()));
    }
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

bool node::contains(kmer const& potential_match) const {
    return subnodes_.find(potential_match.get_hash()) != subnodes_.end();
}

boost::unordered_map<node*, size_t> const& node::get_output_edges() const noexcept {
    return output_edges_;
}

boost::unordered_map<u64, std::vector<subnode>> const& node::get_subnodes() const noexcept {
    return subnodes_;
}

void node::for_every_subnode(std::function<void(subnode const&)> f) const {
    for (auto const& entry : subnodes_) {
        for (subnode const& current : entry.second) {
            f(current);
        }
    }
}

} // namespace pog
