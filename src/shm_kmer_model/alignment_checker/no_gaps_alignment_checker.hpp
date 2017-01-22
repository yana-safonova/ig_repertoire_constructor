//
// Created by Andrew Bzikadze on 5/20/16.
//

#pragma once

#include "shm_config.hpp"
#include "abstract_alignment_checker.hpp"

class NoGapsAlignmentChecker: public ns_abstract_alignment_checker::AbstractAlignmentChecker {
public:
    explicit NoGapsAlignmentChecker(const shm_config::alignment_checker_params &) { }
    virtual bool check(const ns_gene_alignment::EvolutionaryEdgeAlignment &) const;
    virtual ~NoGapsAlignmentChecker() { }
};
