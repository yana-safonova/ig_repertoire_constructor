//
// Created by Andrew Bzikadze on 5/20/16.
//

#pragma once

#include <memory>
#include "evolutionary_edge_alignment/evolutionary_edge_alignment.hpp"

namespace ns_abstract_alignment_checker {
class AbstractAlignmentChecker {
public:
    virtual bool check(const ns_gene_alignment::EvolutionaryEdgeAlignment &) const = 0;
    virtual ~AbstractAlignmentChecker() { }
};
using AbstractAlignmentCheckerPtr = std::shared_ptr<AbstractAlignmentChecker>;
}