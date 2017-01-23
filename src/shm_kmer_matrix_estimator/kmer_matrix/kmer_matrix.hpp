//
// Created by Andrew Bzikadze on 11/11/16.
//

#pragma once

#include <array>

#include "seqan/basic.h"

#include "kmer_utils/kmer_indexed_vector.hpp"

namespace shm_kmer_matrix_estimator {

using KmerMatrix = KmerIndexedVector<std::array<unsigned int, seqan::ValueSize<seqan::Dna>::VALUE>>;

} // End namespace shm_kmer_matrix_estimator