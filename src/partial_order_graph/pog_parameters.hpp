#pragma once

#include <cstring> // size_t
#include "seqan_read.hpp"

namespace pog {

using u64 = unsigned long long;

struct pog_parameters {
    static pog_parameters& instance();

    void set_kmer_size(size_t kmer_size);
    size_t get_kmer_size() const noexcept;

    static u64 const alphabet_size = static_cast<u64>(seqan::ValueSize<nt_t>::VALUE);
    u64 mask;
    float mismatch_penalty;
    float gap_penalty;

private:

    pog_parameters();
    pog_parameters(pog_parameters const&) = delete;
    pog_parameters(pog_parameters&&) = delete;
    size_t kmer_size_;

};

} // namespace pog;
