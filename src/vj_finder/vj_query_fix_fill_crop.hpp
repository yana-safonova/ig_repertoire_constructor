#pragma once

#include "vj_finder_config.hpp"
#include "vj_alignment_structs.hpp"

namespace vj_finder {
    class BaseFillFixCropProcessor {
    protected:
        const VJFinderConfig::AlgorithmParams::FixCropFillParams &params_;
        core::ReadArchive &read_archive_;

    public:
        BaseFillFixCropProcessor(const VJFinderConfig::AlgorithmParams::FixCropFillParams &params,
                                 core::ReadArchive &read_archive) : params_(params),
                                                                    read_archive_(read_archive) { }

        virtual VJHits Process(VJHits vj_hits) = 0;

        virtual ~BaseFillFixCropProcessor() { }
    };

    class AggressiveFillFixCropProcessor : public BaseFillFixCropProcessor {
    public:
        AggressiveFillFixCropProcessor(const VJFinderConfig::AlgorithmParams::FixCropFillParams &params,
                                       core::ReadArchive &read_archive) : BaseFillFixCropProcessor(params,
                                                                                                    read_archive) { }


        VJHits Process(VJHits vj_hits);
    };
}