#pragma once

#include "vj_finder_config.hpp"
#include "vj_alignment_structs.hpp"

namespace vj_finder {
    class FillFixCropProcessor {
        const VJFinderConfig::AlgorithmParams::FixCropFillParams &params_;

        VJHits PerformFixing(VJHits vj_hits);

        VJHits PerformFillingCropping(VJHits vj_hits);

    public:
        FillFixCropProcessor(const VJFinderConfig::AlgorithmParams::FixCropFillParams &params) :
                params_(params) { }

        VJHits Process(VJHits vj_hits);
    };
}