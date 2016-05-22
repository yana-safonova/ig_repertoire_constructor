//
// Created by Andrew Bzikadze on 5/22/16.
//

#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

namespace ns_abstract_mutation_strategy {
using MutationsStatistics = std::unordered_map<std::string, std::vector<unsigned int>>;

class AbstractMutationStrategy {
public:
    virtual MutationsStatistics calculate_mutation_statistics() = 0;
    virtual ~AbstractMutationStrategy() { }
};
using AbstractMutationStrategyPtr = std::shared_ptr<AbstractMutationStrategy>;
}