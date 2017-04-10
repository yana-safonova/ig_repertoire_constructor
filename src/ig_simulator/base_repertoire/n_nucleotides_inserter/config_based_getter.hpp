//
// Created by Andrew Bzikadze on 3/31/17.
//

#pragma once

#include "abstract_n_nucleotides_inserter.hpp"
#include "ig_simulator_config.hpp"

namespace ig_simulator {

AbstractNNucleotidesInserterCPtr get_nucleotides_inserter(const NNucleotidesInserterParams & config);

} // End namespace ig_simulator
