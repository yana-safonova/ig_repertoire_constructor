//
// Created by Andrew Bzikadze on 3/31/17.
//

#pragma once

#include "metaroot_creator.hpp"

namespace ig_simulator {

AbstractMetarootCreatorCPtr get_metarootcreator(const germline_utils::ChainType chain_type,
                                                const MetarootSimulationParams& config,
                                                std::vector<germline_utils::CustomGeneDatabase>& db);

} // End namespace ig_simulator
