#pragma once

#include "vdj_config.hpp"
#include "read_archive.hpp"
#include "vdj_alignments/vdj_hits_storage.hpp"

namespace vdj_labeler {

class VDJLabelerLaunch {
    const VDJLabelerConfig &config_;
    void TestRecombinationCalculator(const core::ReadArchive& reads_archive,
                                     VDJHitsStoragePtr &hits_storage);

public:
    VDJLabelerLaunch(const VDJLabelerConfig &config) :
            config_(config) { }

    void Launch();
};

} // End namespace vdj_labeler