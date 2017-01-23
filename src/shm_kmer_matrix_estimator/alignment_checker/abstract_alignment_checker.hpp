//
// Created by Andrew Bzikadze on 5/20/16.
//

#pragma once

#include <memory>
#include "evolutionary_edge_alignment/evolutionary_edge_alignment.hpp"

namespace shm_kmer_matrix_estimator {

class AbstractAlignmentChecker {
public:
    virtual bool check(EvolutionaryEdgeAlignment &) const = 0;
    virtual ~AbstractAlignmentChecker() { }
};

using AbstractAlignmentCheckerPtr = std::unique_ptr<AbstractAlignmentChecker>;

} // End namespace shm_kmer_matrix_estimator