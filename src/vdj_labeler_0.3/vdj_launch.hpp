#pragma once

#include "vdj_config.hpp"
#include "read_archive.hpp"

namespace vdj_labeler {

class VDJLabelerLaunch {
    const VDJLabelerConfig &config_;

public:
    VDJLabelerLaunch(const VDJLabelerConfig &config) :
            config_(config) { }

    void Launch();
};

} // End namespace vdj_labeler