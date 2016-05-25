#pragma once

#include "cdr_config.hpp"

namespace cdr_labeler {
    class CDRLabelerLaunch {
        const CDRLabelerConfig &config_;

    public:
        CDRLabelerLaunch(const CDRLabelerConfig &config) :
                config_(config) { }

        void Launch();
    };
}