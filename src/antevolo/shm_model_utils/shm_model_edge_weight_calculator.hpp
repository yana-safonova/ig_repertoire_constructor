//
// Created by Andrew Bzikadze on 3/11/17.
//

#pragma once

#include "mutation_strategies/abstract_mutation_strategy.hpp"
#include "shm_model.hpp"
#include "evolutionary_graph_utils/evolutionary_edge/base_evolutionary_edge.hpp"
#include "evolutionary_edge_alignment/evolutionary_edge_alignment.hpp"

namespace antevolo {

class ShmModelEdgeWeightCalculator {
private:
    ShmModel model_;
    shm_kmer_matrix_estimator::AbstractMutationStrategyPtr mutation_strategy_;

private:
    double calculate_weigth_edge_per_position(
        const shm_kmer_matrix_estimator::EvolutionaryEdgeAlignment &src_dst_pair,
        const size_t center_nucl_pos,
        const size_t kmer_len_) const;

    shm_kmer_matrix_estimator::EvolutionaryEdgeAlignment get_prepared_strings(const BaseEvolutionaryEdge& edge) const;

public:
    ShmModelEdgeWeightCalculator() = delete;

    ShmModelEdgeWeightCalculator(ShmModel model,
                                 shm_kmer_matrix_estimator::AbstractMutationStrategyPtr mutation_strategy) :
        model_(std::move(model)), mutation_strategy_(std::move(mutation_strategy)) { }

    double calculate_weigth_edge(const BaseEvolutionaryEdge& edge) const;
};

} // End namespace antevolo