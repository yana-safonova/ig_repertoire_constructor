#pragma once

#include "antevolo_config.hpp"

namespace antevolo {
    class AntEvoloLaunch {
        const AntEvoloConfig& config_;

    public:
        AntEvoloLaunch(const AntEvoloConfig& config) : config_(config) { }

        void Launch();
    };
}