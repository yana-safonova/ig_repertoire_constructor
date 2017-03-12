//
// Created by Andrew Bzikadze on 5/22/16.
//

#pragma once

#include <memory>
#include <unordered_map>
#include <vector>
#include <string>

#include "shm_kmer_matrix_estimator_config.hpp"
#include "../evolutionary_edge_alignment/evolutionary_edge_alignment.hpp"

namespace shm_kmer_matrix_estimator {

class AbstractMutationStrategy {
protected:
    unsigned int kmer_len_;
public:
    explicit AbstractMutationStrategy(const shm_kmer_matrix_estimator_config::mutations_strategy_params &config) :
        kmer_len_(config.kmer_len) { }

    virtual std::vector<size_t> calculate_relevant_positions(EvolutionaryEdgeAlignment &) const = 0;
    virtual ~AbstractMutationStrategy() { }
};
using AbstractMutationStrategyPtr = std::unique_ptr<AbstractMutationStrategy>;

} // End namespace shm_kmer_matrix_estimator