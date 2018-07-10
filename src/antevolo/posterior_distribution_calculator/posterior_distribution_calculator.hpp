//
// Created by Andrew Bzikadze on 3/10/17.
//

#pragma once

#include "shm_model_utils/shm_model.hpp"
#include "kmer_matrix/kmer_matrix.hpp"

namespace antevolo {

class PosteriorDistributionCalculator {
public:
    ShmModel calculate(const ShmModel&,
                       const shm_kmer_matrix_estimator::KmerMatrix& fr_matrix,
                       const shm_kmer_matrix_estimator::KmerMatrix& cdr_matrix) const;
};

} // End namespace antevolo
