#pragma once

#include "vj_query_aligner.hpp"
#include "vj_alignment_structs.hpp"
#include "vj_finder_config.hpp"
#include "vj_query_fix_fill_crop.hpp"

namespace vj_finder {
    class VJQueryProcessor {
        const VJFinderConfig::AlgorithmParams &params_;
        core::ReadArchive &read_archive_;
        const germline_utils::CustomGeneDatabase &v_db_;
        const germline_utils::CustomGeneDatabase &j_db_;

        ProcessedVJHits ComputeFilteringResults(VJHits vj_hits);

        std::shared_ptr<BaseFillFixCropProcessor> GetFillFixCropProcessor();

    public:
        VJQueryProcessor(const VJFinderConfig::AlgorithmParams &params,
                         core::ReadArchive &read_archive,
                         const germline_utils::CustomGeneDatabase &v_db,
                         const germline_utils::CustomGeneDatabase &j_db) : params_(params),
                                                                           read_archive_(read_archive),
                                                                           v_db_(v_db),
                                                                           j_db_(j_db) { }

        ProcessedVJHits Process(const core::Read &read);
    };
}