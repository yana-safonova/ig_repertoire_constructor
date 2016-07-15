#include "pog_parameters.hpp"

namespace pog {

pog_parameters& pog_parameters::instance() {
    static pog_parameters singletone;
    return singletone;
}

pog_parameters::pog_parameters() {
    set_kmer_size(10);
}

void pog_parameters::set_kmer_size(size_t kmer_size) {
    if (kmer_size == kmer_size_)
        return;

    kmer_size_ = kmer_size;
    mask = 1;
    for (size_t i = 0; i < kmer_size; ++i) {
        mask *= seqan::ValueSize<nt_t>::VALUE;
    }
}

size_t pog_parameters::get_kmer_size() const noexcept {
    return kmer_size_;
}

} // namespace pog
