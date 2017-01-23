//
// Created by Andrew Bzikadze on 5/24/16.
//

#pragma once

#include "shm_kmer_matrix_estimator_config.hpp"
#include "abstract_mutation_strategy.hpp"

namespace shm_kmer_matrix_estimator {

class NoKNeighboursMutationStrategy: public AbstractMutationStrategy {
public:
    explicit NoKNeighboursMutationStrategy(const shm_kmer_matrix_estimator_config::mutations_strategy_params &config) :
        AbstractMutationStrategy(config) { }

    std::vector<size_t> calculate_relevant_positions(EvolutionaryEdgeAlignment &) const override;
    virtual ~NoKNeighboursMutationStrategy() { }
};

} // End namespace shm_kmer_matrix_estimator
