//
// Created by Andrew Bzikadze on 5/22/16.
//

#pragma once

#include "../shm_config.hpp"
#include "abstract_mutation_strategy.hpp"

class TrivialMutationStrategy: public ns_abstract_mutation_strategy::AbstractMutationStrategy {
public:
    explicit TrivialMutationStrategy(const shm_config::mutations_strategy_params &config) :
        AbstractMutationStrategy(config) { }

    std::vector<size_t> calculate_relevant_positions(ns_gene_alignment::ReadGermlineAlignment &) const;
    virtual ~TrivialMutationStrategy() { }
};