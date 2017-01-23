//
// Created by Andrew Bzikadze on 5/27/16.
//

#pragma once

#include "shm_kmer_matrix_estimator_config.hpp"

namespace shm_kmer_matrix_estimator {

void parse_command_line_args(shm_kmer_matrix_estimator_config &cfg, int argc, char **argv);

} // End namespace shm_kmer_matrix_estimator