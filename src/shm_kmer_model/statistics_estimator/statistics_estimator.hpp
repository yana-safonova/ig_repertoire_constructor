//
// Created by Andrew Bzikadze on 5/24/16.
//

#pragma once

#include "../shm_config.hpp"
#include "../evolutionary_edge_alignment/evolutionary_edge_alignment.hpp"
#include "../mutation_strategies/abstract_mutation_strategy.hpp"
#include "mutation_statistics.hpp"

class StatisticsEstimator {
private:
    ns_abstract_mutation_strategy::AbstractMutationStrategyPtr mutation_strategy_;
    unsigned int kmer_len_;

public:
    explicit StatisticsEstimator(const shm_config::mutations_strategy_params &config);

    MutationsStatistics calculate_mutation_statistics(ns_gene_alignment::VectorEvolutionaryEdgeAlignments &) const;
};
