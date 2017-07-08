//
// Created by Andrew Bzikadze on 5/20/16.
//

#include <string>

#include "no_gaps_alignment_checker.hpp"

namespace shm_kmer_matrix_estimator {

bool NoGapsAlignmentChecker::check(EvolutionaryEdgeAlignment &evolutionary_edge) const {
    if (evolutionary_edge.IsChecked())
        return evolutionary_edge.CheckIsOk();

    evolutionary_edge.SetChecked();
    bool base_check = AbstractAlignmentChecker::check(evolutionary_edge);

    if (not base_check) {
        evolutionary_edge.SetCheckResult(false);
        return false;
    }

    bool check_result = evolutionary_edge.parent().find_first_of('-') == std::string::npos &&
                        evolutionary_edge.son().   find_first_of('-') == std::string::npos;
    evolutionary_edge.SetCheckResult(check_result);
    return check_result;
}

} // End namespace shm_kmer_matrix_estimator