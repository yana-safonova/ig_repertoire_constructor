//
// Created by Andrew Bzikadze on 3/29/17.
//

#pragma once

#include <memory>
#include <vector>
#include <ostream>

#include "ig_simulator_config.hpp"
#include "germline_utils/germline_db_generator.hpp"
#include "metaroot_cluster/metaroot_cluster.hpp"
#include "metaroot_creator/config_based_getter.hpp"
#include "multiplicity_creator/multiplicity_creator.hpp"
#include "base_repertoire/base_repertoire.hpp"
#include "productivity_checker/productivity_checker.hpp"

namespace ig_simulator {

class BaseRepertoireSimulator {
private:
    AbstractMetarootCreatorCPtr metaroot_creator_p;
    AbstractMultiplicityCreatorPtr multiplicity_creator_p;
    double productive_part;

public:
    BaseRepertoireSimulator(const IgSimulatorConfig::SimulationParams::BaseRepertoireParams& config,
                            const germline_utils::ChainType& chain_type,
                            std::vector<germline_utils::CustomGeneDatabase> &db):
        metaroot_creator_p(get_metarootcreator(chain_type, config.metaroot_simulation_params, db)),
        multiplicity_creator_p(get_multiplicity_creator(config.multiplicity_creator_params)),
        productive_part(config.productive_params.productive_part)
    { }

    BaseRepertoireSimulator() = delete;
    BaseRepertoireSimulator(const BaseRepertoireSimulator&) = delete;
    BaseRepertoireSimulator(BaseRepertoireSimulator&&) = default;
    BaseRepertoireSimulator& operator=(const BaseRepertoireSimulator&) = delete;
    BaseRepertoireSimulator& operator=(BaseRepertoireSimulator&&) = delete;

    BaseRepertoire Simulate(size_t size);
};

} // End namespace ig_simulator
