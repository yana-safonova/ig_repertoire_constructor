//
// Created by Andrew Bzikadze on 5/20/16.
//

#pragma once

#include <memory>
#include "evolutionary_edge_alignment/evolutionary_edge_alignment.hpp"
#include "shm_kmer_matrix_estimator_config.hpp"

namespace shm_kmer_matrix_estimator {

class AbstractAlignmentChecker {
protected:
    const FunctionalityMethod functionality_method;

public:
    explicit AbstractAlignmentChecker(const shm_kmer_matrix_estimator_config::alignment_checker_params& config):
        functionality_method(config.functionality_method)
    { }

    virtual bool check(EvolutionaryEdgeAlignment &alignment) const {
        if (functionality_method == FunctionalityMethod::all)
            return true;
        if (functionality_method == FunctionalityMethod::productive)
            return alignment.Productive();
        return not alignment.Productive();
    };
    virtual ~AbstractAlignmentChecker() { }
};

using AbstractAlignmentCheckerPtr = std::unique_ptr<AbstractAlignmentChecker>;

} // End namespace shm_kmer_matrix_estimator