#pragma once

#include "vj_finder_config.hpp"
#include "vj_alignment_structs.hpp"

namespace vj_finder {
    class FillFixCropProcessor {
        const VJFinderConfig::AlgorithmParams::FixCropFillParams &params_;
        core::ReadArchive &read_archive_;

        VJHits PerformFixing(VJHits vj_hits);

        VJHits PerformFillingCropping(VJHits vj_hits);

    public:
        FillFixCropProcessor(const VJFinderConfig::AlgorithmParams::FixCropFillParams &params,
                             core::ReadArchive &read_archive) :
                params_(params), read_archive_(read_archive) { }

        VJHits Process(VJHits vj_hits);
    };
}