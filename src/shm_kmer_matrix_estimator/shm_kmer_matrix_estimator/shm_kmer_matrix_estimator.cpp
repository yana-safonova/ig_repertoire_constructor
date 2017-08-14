//
// Created by Andrew Bzikadze on 5/24/16.
//

#include "mutation_strategies/trivial_strategy.hpp"
#include "mutation_strategies/no_k_neighbours.hpp"
#include "alignment_checker/no_gaps_alignment_checker.hpp"
#include "shm_kmer_matrix_estimator.hpp"

namespace shm_kmer_matrix_estimator {

ShmKmerMatrixEstimator::ShmKmerMatrixEstimator(const shm_kmer_matrix_estimator_config::mutations_strategy_params &shm_config_ms,
                                         const shm_kmer_matrix_estimator_config::alignment_checker_params &shm_config_ach,
                                         const shm_kmer_matrix_estimator_config::alignment_cropper_params &shm_config_acp) :
    kmer_len_(shm_config_ms.kmer_len)
{
    using MutationStrategyMethod = shm_kmer_matrix_estimator_config::mutations_strategy_params::MutationsStrategyMethod;
    if (shm_config_ms.mutations_strategy_method == MutationStrategyMethod::Trivial)
        mutation_strategy_ = std::unique_ptr<TrivialMutationStrategy>(new TrivialMutationStrategy(shm_config_ms));
    else if (shm_config_ms.mutations_strategy_method == MutationStrategyMethod::NoKNeighbours)
        mutation_strategy_ = std::unique_ptr<NoKNeighboursMutationStrategy>(new NoKNeighboursMutationStrategy(shm_config_ms));

    using AlignmentCheckerMethod = shm_kmer_matrix_estimator_config::alignment_checker_params::AlignmentCheckerMethod;
    if (shm_config_ach.alignment_checker_method == AlignmentCheckerMethod::NoGaps) {
        alignment_checker_ = std::unique_ptr<NoGapsAlignmentChecker>(new NoGapsAlignmentChecker(shm_config_ach));
    }

    using AlignmentCropperMethod = shm_kmer_matrix_estimator_config::alignment_cropper_params::AlignmentCropperMethod ;
    if (shm_config_acp.alignment_cropper_method == AlignmentCropperMethod::UptoLastReliableKMer) {
        alignment_cropper_ =
            std::unique_ptr<UptoLastReliableKmerAlignmentCropper>
            (new UptoLastReliableKmerAlignmentCropper(shm_config_acp.rkmp));
    }
}

void ShmKmerMatrixEstimator::calculate_mutation_statistics_per_position(KmerMatrix &mutations_statistics,
                                                                     const size_t center_nucl_pos,
                                                                     const EvolutionaryEdgeAlignment &alignment) const {
    std::string gene_substring = alignment.parent().substr(center_nucl_pos - kmer_len_ / 2, kmer_len_);
    if ((alignment.son()[center_nucl_pos] == 'N') ||
        (gene_substring.find_first_of('N') != std::string::npos)) {
        return;
    }
    size_t position = seqan::ordValue(static_cast<seqan::Dna>(alignment.son()[center_nucl_pos]));
    mutations_statistics.at(gene_substring).at(position)++;
}


std::pair<KmerMatrix, KmerMatrix>
ShmKmerMatrixEstimator::calculate_mutation_statistics(
        const annotation_utils::AnnotatedCloneSet<annotation_utils::AnnotatedClone>& clone_set) const
{
    VectorEvolutionaryEdgeAlignments alignments;
    for (const auto& clone : clone_set) {
        if (not clone.RegionIsEmpty(annotation_utils::StructuralRegion::CDR1) and
            not clone.RegionIsEmpty(annotation_utils::StructuralRegion::CDR2))
        {
            alignments.emplace_back(clone);
        }
    }
    return calculate_mutation_statistics(alignments);
}


std::pair<KmerMatrix, KmerMatrix>
ShmKmerMatrixEstimator::calculate_mutation_statistics(VectorEvolutionaryEdgeAlignments &alignments) const {
    size_t kmer_matrix_size = static_cast<size_t>(pow(4., kmer_len_));
    KmerMatrix mutations_statistics_fr(kmer_matrix_size);
    KmerMatrix mutations_statistics_cdr(kmer_matrix_size);
    for (auto alignment : alignments) {
        if (not alignment_checker_->check(alignment))
            continue;

        alignment_cropper_->crop(alignment);

        std::vector<size_t> relevant_positions = mutation_strategy_->calculate_relevant_positions(alignment);
        size_t i = 0;

        while (i < relevant_positions.size() and relevant_positions[i] < alignment.cdr1_start()) {
            calculate_mutation_statistics_per_position(mutations_statistics_fr, relevant_positions[i], alignment);
            i++;
        }
        while (i < relevant_positions.size() and relevant_positions[i] <= alignment.cdr1_end()) {
            calculate_mutation_statistics_per_position(mutations_statistics_cdr, relevant_positions[i], alignment);
            i++;
        }
        while (i < relevant_positions.size() and relevant_positions[i] < alignment.cdr2_start()) {
            calculate_mutation_statistics_per_position(mutations_statistics_fr, relevant_positions[i], alignment);
            i++;
        }
        while (i < relevant_positions.size() and relevant_positions[i] <= alignment.cdr2_end()) {
            calculate_mutation_statistics_per_position(mutations_statistics_cdr, relevant_positions[i], alignment);
            i++;
        }
        while (i < relevant_positions.size()) {
            calculate_mutation_statistics_per_position(mutations_statistics_fr, relevant_positions[i], alignment);
            i++;
        }
    }
    return std::make_pair(mutations_statistics_fr, mutations_statistics_cdr);
}

} // End namespace shm_kmer_matrix_estimator
