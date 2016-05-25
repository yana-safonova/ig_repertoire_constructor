//
// Created by Andrew Bzikadze on 5/22/16.
//

#pragma once

#include <memory>
#include <unordered_map>
#include <vector>
#include <string>

#include "../gene_alignment/gene_alignment.hpp"

namespace ns_abstract_mutation_strategy {

class AbstractMutationStrategy {
protected:
    unsigned int kmer_len_;
public:
    explicit AbstractMutationStrategy(const shm_config::mutations_strategy_params &config) :
        kmer_len_(config.kmer_len) { }

    virtual std::vector<size_t> calculate_relevant_positions(ns_gene_alignment::ReadGermlineAlignment &) const = 0;
    virtual ~AbstractMutationStrategy() { }
};
using AbstractMutationStrategyPtr = std::shared_ptr<AbstractMutationStrategy>;
}