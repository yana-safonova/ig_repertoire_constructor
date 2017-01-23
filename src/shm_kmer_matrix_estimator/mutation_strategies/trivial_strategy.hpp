//
// Created by Andrew Bzikadze on 5/22/16.
//

#pragma once

#include "shm_kmer_matrix_estimator_config.hpp"
#include "abstract_mutation_strategy.hpp"

namespace shm_kmer_matrix_estimator {

class TrivialMutationStrategy: public AbstractMutationStrategy {
public:
    explicit TrivialMutationStrategy(const shm_kmer_matrix_estimator_config::mutations_strategy_params &config) :
        AbstractMutationStrategy(config) {}

    virtual std::vector<size_t> calculate_relevant_positions(EvolutionaryEdgeAlignment &) const override;
    virtual ~TrivialMutationStrategy() {}
};

} // End namespace shm_kmer_matrix_estimator