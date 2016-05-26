//
// Created by Andrew Bzikadze on 5/24/16.
//

#pragma once

#include "../shm_config.hpp"
#include "abstract_mutation_strategy.hpp"

class NoKNeighboursMutationStrategy: public ns_abstract_mutation_strategy::AbstractMutationStrategy {
public:
    explicit NoKNeighboursMutationStrategy(const shm_config::mutations_strategy_params &config)
        : AbstractMutationStrategy(config) { }

    std::vector<size_t> calculate_relevant_positions(ns_gene_alignment::ReadGermlineAlignment &) const;
    virtual ~NoKNeighboursMutationStrategy() { }
};
