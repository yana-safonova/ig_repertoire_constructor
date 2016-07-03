#include "vj_hits_filter.hpp"

namespace vj_finder {
    void CustomVjHitsFilter::Add(std::shared_ptr<VJHitsFilter> vj_filter_ptr) {
        vj_filters_.push_back(vj_filter_ptr);
    }

    bool CustomVjHitsFilter::Filter(const VJHits &vj_hits) const {
        for(auto it = vj_filters_.begin(); it != vj_filters_.end(); it++)
            if((*it)->Filter(vj_hits))
                return true;
        return false;
    }

    void VersatileVjFilter::InitializeCustomVjFinder() {
        custom_vj_filter_.Add(std::shared_ptr<VJHitsFilter>(new EmptyHitsFilter()));
        custom_vj_filter_.Add(std::shared_ptr<VJHitsFilter>(new LeftCoverageFilter(
                filtering_params_.left_uncovered_limit)));
        custom_vj_filter_.Add(std::shared_ptr<VJHitsFilter>(new RightCoverageFilter(
                filtering_params_.right_uncovered_limit)));
        custom_vj_filter_.Add(std::shared_ptr<VJHitsFilter>(new VSegmentLengthFilter(
                filtering_params_.min_v_segment_length)));
        custom_vj_filter_.Add(std::shared_ptr<VJHitsFilter>(new JSegmentLengthFilter(
                filtering_params_.min_j_segment_length)));
        custom_vj_filter_.Add(std::shared_ptr<VJHitsFilter>(new AlignedSegmentLengthFilter(
                filtering_params_.min_aligned_length)));
    }
}