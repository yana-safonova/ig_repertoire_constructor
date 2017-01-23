//
// Created by Andrew Bzikadze on 5/24/16.
//

#pragma once

#include "shm_config.hpp"
#include "evolutionary_edge_alignment/evolutionary_edge_alignment.hpp"
#include "mutation_strategies/abstract_mutation_strategy.hpp"
#include "kmer_utils/kmer_indexed_vector.hpp"
#include "kmer_matrix/kmer_matrix.hpp"

namespace shm_kmer_matrix_estimator {

class StatisticsEstimator {
private:
    AbstractMutationStrategyPtr mutation_strategy_;
    unsigned int kmer_len_;

private:
    void calculate_mutation_statistics_per_position(KmerMatrix &mutations_statistics,
                                                    const size_t center_nucl_pos,
                                                    const EvolutionaryEdgeAlignment &) const;

public:
    explicit StatisticsEstimator(const shm_config::mutations_strategy_params &config);


    std::pair<KmerMatrix, KmerMatrix>
    calculate_mutation_statistics(VectorEvolutionaryEdgeAlignments &) const;
};

} // End namespace shm_kmer_matrix_estimator
