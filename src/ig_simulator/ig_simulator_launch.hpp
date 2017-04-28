//
// Created by Andrew Bzikadze on 3/15/17.
//

#pragma once

#include "ig_simulator_config.hpp"
#include "germline_utils/chain_type.hpp"
#include "base_repertoire/base_repertoire.hpp"
#include "clonal_trees/forest/forest.hpp"

namespace ig_simulator {

class IgSimulatorLaunch {
private:
    IgSimulatorConfig config_;

private:
    germline_utils::ChainType GetLaunchChainType() const;

    std::vector<germline_utils::CustomGeneDatabase>
    GetDB(const germline_utils::ChainType chain_type) const;

    BaseRepertoire
    GetBaseRepertoire(const germline_utils::ChainType chain_type,
                      std::vector<germline_utils::CustomGeneDatabase>& db) const;

    template<class PoolManager>
    ForestStorage __GetForestStorage(const BaseRepertoire& base_repertoire) const;

    ForestStorage GetForestStorage(const BaseRepertoire& base_repertoire) const;

public:
    IgSimulatorLaunch(const IgSimulatorConfig &config) :
        config_(config)
    { }

    void Run();

    IgSimulatorLaunch() = delete;
    IgSimulatorLaunch(const IgSimulatorLaunch&) = delete;
    IgSimulatorLaunch(IgSimulatorLaunch&&) = delete;
    IgSimulatorLaunch& operator=(const IgSimulatorLaunch&) = delete;
    IgSimulatorLaunch& operator=(IgSimulatorLaunch&&) = delete;
};

} // End namespace ig_simulator
