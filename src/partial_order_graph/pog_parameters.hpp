#pragma once

#include <cstring> // size_t
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <logger/logger.hpp>
#include <logger/log_writers.hpp>

namespace pog {

using u64 = unsigned long long;

using id_t = seqan::CharString;
using nt_t = seqan::Dna5;
using seq_t = seqan::String<nt_t>;

struct pog_parameters {
    static pog_parameters& instance();

    void set_kmer_size(size_t kmer_size);
    size_t get_kmer_size() const noexcept;

    static u64 const alphabet_size = static_cast<u64>(seqan::ValueSize<nt_t>::VALUE);
    u64 mask;
    float mismatch_penalty = -1.4f;
    float gap_penalty = -2.f;

    float bulge_coverage_difference = 2.f;

private:

    pog_parameters();
    pog_parameters(pog_parameters const&) = delete;
    pog_parameters(pog_parameters&&) = delete;
    size_t kmer_size_;

};

} // namespace pog;
