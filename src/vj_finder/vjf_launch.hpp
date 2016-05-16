#pragma once

#include "vj_finder_config.hpp"

namespace vj_finder {
    class VJFinderLaunch {
        const vjf_config &config_;



    public:
        VJFinderLaunch(const vjf_config &config) :
                config_(config) { }

        void Run();

    private:
        DECL_LOGGER("VJFinderLaunch");
    };
}