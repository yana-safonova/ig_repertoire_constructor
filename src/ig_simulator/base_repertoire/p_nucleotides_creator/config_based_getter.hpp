//
// Created by Andrew Bzikadze on 3/31/17.
//

#pragma once

#include "abstract_nucleotides_creator.hpp"
#include "ig_simulator_config.hpp"

namespace ig_simulator {

AbstractPNucleotidesCreatorCPtr get_nucleotides_creator(const PNucleotidesCreatorParams &config);

} // End namespace ig_simulator
