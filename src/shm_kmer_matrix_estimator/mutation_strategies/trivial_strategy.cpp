//
// Created by Andrew Bzikadze on 5/23/16.
//

#include <string>
#include <seqan/file.h>

#include "trivial_strategy.hpp"

namespace shm_kmer_matrix_estimator {

std::vector<size_t>
TrivialMutationStrategy::calculate_relevant_positions(EvolutionaryEdgeAlignment &alignment) const {
    std::vector<size_t> relevant_positions(alignment.son().size() - kmer_len_ + 1);
    std::iota(relevant_positions.begin(), relevant_positions.end(), kmer_len_ / 2);
    return relevant_positions;
}

} // End namespace shm_kmer_matrix_estimator

