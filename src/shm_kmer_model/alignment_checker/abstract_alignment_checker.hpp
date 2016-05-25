//
// Created by Andrew Bzikadze on 5/20/16.
//

#pragma once

#include <memory>
#include "gene_alignment/gene_alignment.hpp"

namespace ns_abstract_alignment_checker {
class AbstractAlignmentChecker {
public:
    virtual bool check(const ns_gene_alignment::ReadGermlineAlignment&) const = 0;
    virtual ~AbstractAlignmentChecker() { }
};
using AbstractAlignmentCheckerPtr = std::shared_ptr<AbstractAlignmentChecker>;
}