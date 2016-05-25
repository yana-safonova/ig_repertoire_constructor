//
// Created by Andrew Bzikadze on 5/20/16.
//

#include <string>

#include "no_gaps_alignment_checker.hpp"

bool NoGapsAlignmentChecker::check(const ns_gene_alignment::ReadGermlineAlignment &germline_read_pair) const {
    return germline_read_pair.read().find_first_of('-') == std::string::npos &&
        germline_read_pair.germline().find_first_of('-') == std::string::npos;
}