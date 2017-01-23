//
// Created by Andrew Bzikadze on 5/25/16.
//

#include "mutation_statistics.hpp"

namespace shm_kmer_matrix_estimator {

MutationsStatistics::MutationsStatistics(unsigned int kmer_len) {
    std::vector<std::string> kmers;
    generate_kmer_keys(kmers, kmer_len);

    auto alphabet_size = seqan::ValueSize<seqan::Dna>::VALUE;
    for (auto it = kmers.begin(); it != kmers.end(); ++it)
        statistics[*it] = std::vector<unsigned int>(alphabet_size);
}

void MutationsStatistics::generate_kmer_keys(std::vector<std::string> &kmers, const unsigned int kmer_len) const {
    std::string kmer(kmer_len, ' ');
    generate_kmer_keys_(kmers, kmer_len, 0, kmer);
}

void MutationsStatistics::generate_kmer_keys_(std::vector<std::string> &kmers,
                                              const unsigned int kmer_len,
                                              unsigned int curr_len,
                                              std::string kmer) const {
    if (curr_len == kmer_len) {
        kmers.emplace_back(kmer);
        return;
    }
    auto alphabet_size = seqan::ValueSize<seqan::Dna>::VALUE;
    for (size_t i = 0; i < alphabet_size; ++i) {
        kmer[curr_len] = seqan::Dna(i);
        generate_kmer_keys_(kmers, kmer_len, curr_len + 1, kmer);
    }
}

} // End namespace shm_kmer_matrix_estimator
