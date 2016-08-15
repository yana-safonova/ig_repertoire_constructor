#include "vj_query_processing.hpp"
#include "vj_hits_filter.hpp"
#include "vj_query_fix_fill_crop.hpp"

namespace vj_finder {
    std::shared_ptr<BaseFillFixCropProcessor> VJQueryProcessor::GetFillFixCropProcessor() {
        return std::shared_ptr<BaseFillFixCropProcessor>(
                new AggressiveFillFixCropProcessor(params_.fix_crop_fill_params,
                                                   read_archive_));
    }

    ProcessedVJHits VJQueryProcessor::ComputeFilteringResults(const core::Read &read, VJHits vj_hits) {
        ProcessedVJHits processed_hits(read);
        processed_hits.vj_hits = vj_hits;
        if(params_.filtering_params.enable_filtering) {
            VersatileVjFilter vj_filter(params_.filtering_params);
            processed_hits.filtering_info = vj_filter.Filter(vj_hits);
        }
        return processed_hits;
    }

    ProcessedVJHits VJQueryProcessor::Process(const core::Read &read) {
        VJQueryAligner vj_query_aligner(params_, read_archive_, v_db_, j_db_);
        VJHits vj_hits = vj_query_aligner.Align(read);
        ProcessedVJHits hits_after_fitering = ComputeFilteringResults(read, vj_hits);
        if(hits_after_fitering.ReadToBeFiltered()) {
            return hits_after_fitering;
        }
        auto fix_fill_crop_processor = GetFillFixCropProcessor();
        hits_after_fitering.vj_hits = fix_fill_crop_processor->Process(hits_after_fitering.vj_hits);
        return hits_after_fitering;
    }
}