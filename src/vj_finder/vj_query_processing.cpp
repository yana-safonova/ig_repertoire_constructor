#include "vj_query_processing.hpp"
#include "vj_hits_filter.hpp"
#include "vj_query_fix_fill_crop.hpp"

namespace vj_finder {
    ProcessedVJHits VJQueryProcessor::ComputeFilteringResults(VJHits vj_hits) {
        bool read_to_be_filtered = false;
        if(params_.filtering_params.enable_filtering) {
            VersatileVjFilter vj_filter(params_.filtering_params);
            if(vj_filter.Filter(vj_hits))
                read_to_be_filtered = true;
        }
        if(read_to_be_filtered) {
            return ProcessedVJHits();
        }
        return ProcessedVJHits(vj_hits);
    }

    ProcessedVJHits VJQueryProcessor::Process(const core::Read &read) {
        VJQueryAligner vj_query_aligner(params_, v_db_, j_db_);
        VJHits vj_hits = vj_query_aligner.Align(read);
        ProcessedVJHits hits_after_fitering = ComputeFilteringResults(vj_hits);
        if(!hits_after_fitering) {
            return hits_after_fitering;
        }
        FillFixCropProcessor fix_processor(params_.fix_crop_fill_params);
        return ProcessedVJHits(fix_processor.Process(*hits_after_fitering));
    }
}