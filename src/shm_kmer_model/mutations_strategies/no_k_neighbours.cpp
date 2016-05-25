//
// Created by Andrew Bzikadze on 5/24/16.
//

#include <deque>
#include "no_k_neighbours.hpp"

std::vector<size_t> NoKNeighboursMutationStrategy::calculate_relevant_positions
    (ns_gene_alignment::ReadGermlineAlignment &alignment) const
{
    std::deque<int> mismatch_positions(alignment.read().size());
    for (size_t i = 0; i < alignment.read().size(); ++i)
        mismatch_positions[i] = static_cast<int>(alignment.read()[i] != alignment.germline()[i]);

    for (size_t i = 0; i < kmer_len_ / 2; ++i) {
        mismatch_positions.push_back(0);
        mismatch_positions.push_front(0);
    }

    std::vector<int> cum_sums(alignment.read().size());
    for (size_t i = 0; i <= kmer_len_; ++i)
        cum_sums[0] += mismatch_positions[i];

    for (size_t i = 1; i < alignment.read().size(); ++i)
        cum_sums[i] = cum_sums[i - 1] - mismatch_positions[i - 1] + mismatch_positions[i + kmer_len_ - 1];

    std::vector<size_t> relevant_positions;
    for (size_t i = 0; i < alignment.read().size(); ++i) {
        if (cum_sums[i] == 0 ||
            ((cum_sums[i] == 1) && (alignment.read()[i] != alignment.germline()[i])))
            relevant_positions.push_back(i);
    }

    return relevant_positions;
}
