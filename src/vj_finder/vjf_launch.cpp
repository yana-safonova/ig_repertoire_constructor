#include <logger/logger.hpp>

#include "vjf_launch.hpp"

namespace vj_finder {
    void VJFinderLaunch::Run() {
        INFO("== VJ Finder starts == ");
        core::ReadArchive read_archive(config_.iop.input.input_reads);
        INFO("== VJ Finder ends == ");
    }
}