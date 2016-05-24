//
// Created by Andrew Bzikadze on 5/24/16.
//

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace ns_mutation_statistics {
    using MutationsStatistics = std::unordered_map<std::string, std::vector<unsigned int>>;
}