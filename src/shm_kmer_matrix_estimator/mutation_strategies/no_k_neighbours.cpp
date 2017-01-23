//
// Created by Andrew Bzikadze on 5/24/16.
//

#include <deque>
#include "no_k_neighbours.hpp"

namespace shm_kmer_matrix_estimator {

std::vector<size_t>
NoKNeighboursMutationStrategy::calculate_relevant_positions(EvolutionaryEdgeAlignment &alignment) const
{
    std::deque<int> mismatch_positions(alignment.son().size());
    for (size_t i = 0; i < alignment.size(); ++i)
        mismatch_positions[i] = static_cast<int>(alignment.son()[i] != alignment.parent()[i]);

    for (size_t i = 0; i < kmer_len_ / 2; ++i) {
        mismatch_positions.push_back(0);
        mismatch_positions.push_front(0);
    }

    std::vector<int> cum_sums(alignment.son().size());
    for (size_t i = 0; i <= kmer_len_; ++i)
        cum_sums[0] += mismatch_positions[i];

    for (size_t i = 1; i < alignment.size(); ++i)
        cum_sums[i] = cum_sums[i - 1] - mismatch_positions[i - 1] + mismatch_positions[i + kmer_len_ - 1];

    std::vector<size_t> relevant_positions;
    // First and last kmer_len/2 nucleotides are not taken.
    for (size_t i = kmer_len_ / 2; i < alignment.size() - kmer_len_ / 2; ++i) {
        // We take the position `i', if
        //     1) there is no mutation ( == 0);
        // OR
        //     2) there is mutation AND mutation is in the central nucleotide.
        if (cum_sums[i] == 0 ||
            ((cum_sums[i] == 1) && (alignment.son()[i] != alignment.parent()[i])))
            relevant_positions.push_back(i);
    }

    return relevant_positions;
}

} // End namespace shm_kmer_matrix_estimator
