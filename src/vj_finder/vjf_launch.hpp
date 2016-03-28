#include <logger/logger.hpp>

#include "vj_finder_config.hpp"

#include <read_archive.hpp>

namespace vj_finder {
    class VJFinderLaunch {
        const vjf_config &config_;
    public:
        VJFinderLaunch(const vjf_config &config) :
                config_(config) { }

        void Run();
    };
}