//
// Created by Andrew Bzikadze on 5/24/16.
//

#include "../mutation_strategies/trivial_strategy.hpp"
#include "../mutation_strategies/no_k_neighbours.hpp"
#include "statistics_estimator.hpp"

using namespace ns_gene_alignment;

StatisticsEstimator::StatisticsEstimator(const shm_config::mutations_strategy_params &config) :
    kmer_len_(config.kmer_len) {
    using MutationStrategyMethod = shm_config::mutations_strategy_params::MutationsStrategyMethod;
    if (config.mutations_strategy_method == MutationStrategyMethod::Trivial)
        mutation_strategy_ = std::make_shared<TrivialMutationStrategy>(TrivialMutationStrategy(config));
    else if (config.mutations_strategy_method == MutationStrategyMethod::NoKNeighbours)
        mutation_strategy_ = std::make_shared<NoKNeighboursMutationStrategy>(NoKNeighboursMutationStrategy(config));
}

MutationsStatistics StatisticsEstimator::calculate_mutation_statistics(VectorReadGermlineAlignments &alignments) const {
    MutationsStatistics mutations_statistics(kmer_len_);
    for (auto &alignment : alignments) {
        std::vector<size_t> relevant_positions = mutation_strategy_->calculate_relevant_positions(alignment);

        for (auto it = relevant_positions.begin(); it != relevant_positions.end(); ++it) {
            size_t &current_pos = *it;
            size_t center_nucl_pos = current_pos + kmer_len_ / 2;
            std::string gene_substring = alignment.germline().substr(current_pos, kmer_len_);

            if ((alignment.read()[center_nucl_pos] == 'N') ||
                (gene_substring.find_first_of('N') != std::string::npos))
                continue;

            size_t position = seqan::ordValue(static_cast<seqan::Dna>(alignment.read()[center_nucl_pos]));
            mutations_statistics.at(gene_substring).at(position)++;
        }
    }
    return mutations_statistics;
}
