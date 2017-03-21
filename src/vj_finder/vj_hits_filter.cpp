#include "vj_hits_filter.hpp"

namespace vj_finder {
    std::ostream& operator<<(std::ostream& out, const VJFilteringReason &filtering_reason) {
        if(filtering_reason == VJFilteringReason::LeftUncoveredLimitReason) {
            out << "LEFT_UNCOVERED_LIMIT_EXCEEDS";
            return out;
        }
        if(filtering_reason == VJFilteringReason::RightUncoveredLimitReason) {
            out << "RIGHT_UNCOVERED_LIMIT_EXCEEDS";
            return out;
        }
        if(filtering_reason == VJFilteringReason::VJHitsAreEmptyReason) {
            out << "VJ_HITS_ARE_EMPTY";
            return out;
        }
        if(filtering_reason == VJFilteringReason::VSegmentLengthReason) {
            out << "SHORT_V_ALIGNMENT";
            return out;
        }
        if(filtering_reason == VJFilteringReason::JSegmentLengthReason) {
            out << "SHORT_J_ALIGNMENT";
            return out;
        }
        if(filtering_reason == VJFilteringReason::AlignedSegmentLengthReason) {
            out << "SHORT_ALIGNED_SEGMENT";
            return out;
        }
        VERIFY_MSG(false, "Unknown filtering reason");
        return out;
    }

    std::ostream& operator<<(std::ostream& out, const VJFilteringInfo &filtering_info) {
        if(filtering_info.filtering_reason == VJFilteringReason::LeftUncoveredLimitReason) {
            out << filtering_info.filtering_reason << "\t" << filtering_info.left_uncovered_length;
            return out;
        }
        if(filtering_info.filtering_reason == VJFilteringReason::RightUncoveredLimitReason) {
            out << filtering_info.filtering_reason << "\t" << filtering_info.right_uncovered_length;
            return out;
        }
        if(filtering_info.filtering_reason == VJFilteringReason::VJHitsAreEmptyReason) {
            out << filtering_info.filtering_reason;
            return out;
        }
        if(filtering_info.filtering_reason == VJFilteringReason::VSegmentLengthReason) {
            out << filtering_info.filtering_reason << "\t" << filtering_info.v_segment_length;
            return out;
        }
        if(filtering_info.filtering_reason == VJFilteringReason::JSegmentLengthReason) {
            out << filtering_info.filtering_reason << "\t" << filtering_info.j_segment_length;
            return out;
        }
        if(filtering_info.filtering_reason == VJFilteringReason::AlignedSegmentLengthReason) {
            out << filtering_info.filtering_reason << "\t" << filtering_info.aligned_segment_length;
            return out;
        }
        VERIFY_MSG(false, "Unknown filtering reason");
        return out;
    }

    void CustomVjHitsFilter::Add(std::shared_ptr<VJHitsFilter> vj_filter_ptr) {
        vj_filters_.push_back(vj_filter_ptr);
    }

    VJFilteringInfo CustomVjHitsFilter::Filter(const VJHits &vj_hits) const {
        for(auto it = vj_filters_.begin(); it != vj_filters_.end(); it++) {
            auto filtering_res = (*it)->Filter(vj_hits);
            if(filtering_res.read_to_be_filtered)
                return filtering_res;
        }
        return VJFilteringInfo(vj_hits.Read());
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