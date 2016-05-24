//
// Created by Andrew Bzikadze on 5/20/16.
//

#include <string>

#include "no_gaps_alignment_checker.hpp"
#include "gene_alignment/gene_alignment.hpp"

bool NoGapsAlignmentChecker::check(const ns_gene_alignment::ReadGermlinePair & germline_read_pair) {
    return germline_read_pair.first.find_first_of('-') == std::string::npos &&
           germline_read_pair.second.find_first_of('-') == std::string::npos;
}