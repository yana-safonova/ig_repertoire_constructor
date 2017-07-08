//
// Created by Andrew Bzikadze on 5/20/16.
//

#pragma once

#include "shm_kmer_matrix_estimator_config.hpp"
#include "abstract_alignment_checker.hpp"

namespace shm_kmer_matrix_estimator {

class NoGapsAlignmentChecker: public AbstractAlignmentChecker {
public:
    explicit NoGapsAlignmentChecker(const shm_kmer_matrix_estimator_config::alignment_checker_params &config) :
        AbstractAlignmentChecker(config)
    {}

    virtual bool check(EvolutionaryEdgeAlignment &) const override;
    virtual ~NoGapsAlignmentChecker() { }
};

} // End namespace shm_kmer_matrix_estimator
