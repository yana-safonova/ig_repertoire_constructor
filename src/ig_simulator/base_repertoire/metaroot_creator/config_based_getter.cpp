//
// Created by Andrew Bzikadze on 3/31/17.
//

#include "config_based_getter.hpp"


namespace ig_simulator {

AbstractMetarootCreatorCPtr get_metarootcreator(const germline_utils::ChainType chain_type,
                                               const MetarootSimulationParams& config,
                                               std::vector<germline_utils::CustomGeneDatabase>& db)
{
    if (chain_type.IsVDJ())
        return AbstractMetarootCreatorCPtr(new VDJMetarootCreator(config, db));
    return AbstractMetarootCreatorCPtr(new VJMetarootCreator(config, db));
}

} // End namespace ig_simulator

