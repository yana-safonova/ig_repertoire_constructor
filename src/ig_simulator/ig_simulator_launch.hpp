//
// Created by Andrew Bzikadze on 3/15/17.
//

#pragma once

#include "ig_simulator_config.hpp"
#include "germline_utils/chain_type.hpp"

namespace ig_simulator {

class IgSimulatorLaunch {
private:
    IgSimulatorConfig config_;

private:
    germline_utils::ChainType GetLaunchChainType() const;

public:
    IgSimulatorLaunch(const IgSimulatorConfig &config) :
        config_(config)
    { }

    IgSimulatorLaunch(const IgSimulatorLaunch&) = delete;
    IgSimulatorLaunch(IgSimulatorLaunch&&) = delete;

    void Run();
};

} // End namespace ig_simulator
