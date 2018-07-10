//
// Created by Andrew Bzikadze on 5/24/16.
//

#pragma once

#include "shm_kmer_matrix_estimator_config.hpp"
#include "evolutionary_edge_alignment/evolutionary_edge_alignment.hpp"
#include "mutation_strategies/abstract_mutation_strategy.hpp"
#include "alignment_checker/abstract_alignment_checker.hpp"
#include "alignment_cropper/abstract_alignment_cropper.hpp"
#include "alignment_cropper/upto_last_reliable_kmer_alignment_cropper.hpp"
#include "kmer_utils/kmer_indexed_vector.hpp"
#include "kmer_matrix/kmer_matrix.hpp"
#include "annotation_utils/annotated_clone_set.hpp"

namespace shm_kmer_matrix_estimator {

class ShmKmerMatrixEstimator {
private:
    AbstractMutationStrategyPtr mutation_strategy_;
    AbstractAlignmentCheckerPtr alignment_checker_;
    AbstractAlignmentCropperPtr alignment_cropper_;
    unsigned int kmer_len_;

private:
    void calculate_mutation_statistics_per_position(KmerMatrix &mutations_statistics,
                                                    const size_t center_nucl_pos,
                                                    const EvolutionaryEdgeAlignment &) const;

public:
    explicit ShmKmerMatrixEstimator(const shm_kmer_matrix_estimator_config::mutations_strategy_params &shm_config_ms,
                                    const shm_kmer_matrix_estimator_config::alignment_checker_params &shm_config_ach,
                                    const shm_kmer_matrix_estimator_config::alignment_cropper_params &shm_config_acp);


    std::pair<KmerMatrix, KmerMatrix>
    calculate_mutation_statistics(VectorEvolutionaryEdgeAlignments &) const;

    std::pair<KmerMatrix, KmerMatrix>
    calculate_mutation_statistics(
        const annotation_utils::AnnotatedCloneSet<annotation_utils::AnnotatedClone>& clone_set) const;
};

} // End namespace shm_kmer_matrix_estimator
