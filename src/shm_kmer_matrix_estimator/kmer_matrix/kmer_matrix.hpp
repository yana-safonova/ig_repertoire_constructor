//
// Created by Andrew Bzikadze on 11/11/16.
//

#pragma once

#include <array>

#include "seqan/basic.h"

#include "kmer_utils/kmer_indexed_vector.hpp"

namespace shm_kmer_matrix_estimator {

using KmerMatrixRowType = std::array<unsigned int, seqan::ValueSize<seqan::Dna>::VALUE>;
using KmerMatrix = KmerIndexedVector<KmerMatrixRowType>;

} // End namespace shm_kmer_matrix_estimator