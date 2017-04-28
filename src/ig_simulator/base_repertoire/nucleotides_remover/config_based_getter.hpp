//
// Created by Andrew Bzikadze on 3/31/17.
//

#pragma once

#include "abstract_nucleotides_remover.hpp"
#include "ig_simulator_config.hpp"

namespace ig_simulator {

AbstractNucleotidesRemoverCPtr get_nucleotides_remover(const NucleotidesRemoverParams & config);

} // End namespace ig_simulator
