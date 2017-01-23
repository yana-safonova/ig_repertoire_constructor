//
// Created by Andrew Bzikadze on 5/20/16.
//

#include <string>

#include "no_gaps_alignment_checker.hpp"

namespace shm_kmer_matrix_estimator {

bool NoGapsAlignmentChecker::check(const EvolutionaryEdgeAlignment &germline_read_pair) const {
    return germline_read_pair.parent().find_first_of('-') == std::string::npos &&
           germline_read_pair.son().   find_first_of('-') == std::string::npos;
}

} // End namespace shm_kmer_matrix_estimator