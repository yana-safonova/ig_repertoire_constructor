//
// Created by Andrew Bzikadze on 5/24/16.
//

#include "mutation_strategies/trivial_strategy.hpp"
#include "mutation_strategies/no_k_neighbours.hpp"
#include "alignment_checker/no_gaps_alignment_checker.hpp"
#include "statistics_estimator.hpp"

namespace shm_kmer_matrix_estimator {

StatisticsEstimator::StatisticsEstimator(const shm_config::mutations_strategy_params &config_ms,
                                         const shm_config::alignment_checker_params &config_ach) :
    kmer_len_(config_ms.kmer_len)
{
    using MutationStrategyMethod = shm_config::mutations_strategy_params::MutationsStrategyMethod;
    if (config_ms.mutations_strategy_method == MutationStrategyMethod::Trivial)
        mutation_strategy_ = std::unique_ptr<TrivialMutationStrategy>(new TrivialMutationStrategy(config_ms));
    else if (config_ms.mutations_strategy_method == MutationStrategyMethod::NoKNeighbours)
        mutation_strategy_ = std::unique_ptr<NoKNeighboursMutationStrategy>(new NoKNeighboursMutationStrategy(config_ms));

    using AlignmentCheckerMethod = shm_config::alignment_checker_params::AlignmentCheckerMethod;
    if (config_ach.alignment_checker_method == AlignmentCheckerMethod::NoGaps) {
        alignment_checker_ = std::unique_ptr<NoGapsAlignmentChecker>(new NoGapsAlignmentChecker(config_ach));
    }
}

void StatisticsEstimator::calculate_mutation_statistics_per_position(KmerMatrix &mutations_statistics,
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
StatisticsEstimator::calculate_mutation_statistics(VectorEvolutionaryEdgeAlignments &alignments) const {
    size_t kmer_matrix_size = static_cast<size_t>(pow(4., kmer_len_));
    KmerMatrix mutations_statistics_fr(kmer_matrix_size);
    KmerMatrix mutations_statistics_cdr(kmer_matrix_size);
    for (auto &alignment : alignments) {
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
        while (i < relevant_positions.size() and relevant_positions[i] < alignment.cdr2_end()) {
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
