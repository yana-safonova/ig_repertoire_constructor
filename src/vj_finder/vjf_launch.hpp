#pragma once

#include "vj_finder_config.hpp"

namespace vj_finder {
    class VJFinderLaunch {
        const VJFinderConfig &config_;



    public:
        VJFinderLaunch(const VJFinderConfig &config) :
                config_(config) { }

        void Run();

    private:
        DECL_LOGGER("VJFinderLaunch");
    };
}