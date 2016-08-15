#pragma once

#include "vj_finder_config.hpp"
#include "vj_alignment_structs.hpp"

namespace vj_finder {
    enum VJFilteringReason { UnknownFilteringReason, VJHitsAreEmptyReason, LeftUncoveredLimitReason,
        RightUncoveredLimitReason, VSegmentLengthReason, JSegmentLengthReason, AlignedSegmentLengthReason };

    std::ostream& operator<<(std::ostream& out, const VJFilteringReason &filtering_reason);

    struct VJFilteringInfo {
        const core::Read *read;
        bool read_to_be_filtered; // true if read should be filtered out
        VJFilteringReason filtering_reason;
        int left_uncovered_length;
        int right_uncovered_length;
        size_t v_segment_length;
        size_t j_segment_length;
        size_t aligned_segment_length;

        VJFilteringInfo(const core::Read &read) : read(&read),
                                                    read_to_be_filtered(false),
                                                    filtering_reason(VJFilteringReason::UnknownFilteringReason),
                                                    left_uncovered_length(std::numeric_limits<int>::max()),
                                                    right_uncovered_length(std::numeric_limits<int>::max()),
                                                    v_segment_length(size_t(-1)),
                                                    j_segment_length(size_t(-1)),
                                                    aligned_segment_length(size_t(-1)) { }

        void UpdateFilteringReason(VJFilteringReason filtering_reason) {
            if(filtering_reason != VJFilteringReason::UnknownFilteringReason) {
                read_to_be_filtered = true;
            }
            this->filtering_reason = filtering_reason;
        }
    };

    std::ostream& operator<<(std::ostream& out, const VJFilteringInfo &filtering_info);

    class VJHitsFilter {
    public:
        // returns true if VJ Hits are bad, otherwise returns false
        virtual VJFilteringInfo Filter(const VJHits &vj_hits) const = 0;

        ~VJHitsFilter() { }
    };

    class EmptyHitsFilter : public VJHitsFilter {
    public:
        VJFilteringInfo Filter(const VJHits &vj_hits) const {
            bool empty = vj_hits.NumJHits() == 0 or vj_hits.NumVHits() == 0;
            TRACE("V hits are empty and J hits are empty: " << empty);
            VJFilteringInfo res(vj_hits.Read());
            if(empty) {
                res.UpdateFilteringReason(VJFilteringReason::VJHitsAreEmptyReason);
            }
            return res;
        }

    private:
        DECL_LOGGER("EmptyHitsFilter");
    };

    class LeftCoverageFilter : public VJHitsFilter {
        int left_uncovered_limit_;

    public:
        LeftCoverageFilter(int left_uncovered_limit) :
                left_uncovered_limit_(left_uncovered_limit) { }

        VJFilteringInfo Filter(const VJHits &vj_hits) const {
            bool bad = vj_hits.GetVHitByIndex(0).LeftUncovered() > left_uncovered_limit_;
            TRACE("Left uncovered (" << vj_hits.GetVHitByIndex(0).LeftUncovered() <<
                    ") > limit (" << left_uncovered_limit_ << "): " << bad);
            VJFilteringInfo res(vj_hits.Read());
            if(bad) {
                res.UpdateFilteringReason(VJFilteringReason::LeftUncoveredLimitReason);
                res.left_uncovered_length = vj_hits.GetVHitByIndex(0).LeftUncovered();
            }
            return res;
        }

    private:
        DECL_LOGGER("LeftCoverageFilter");
    };

    class RightCoverageFilter : public VJHitsFilter {
        int right_uncovered_limit_;

    public:
        RightCoverageFilter(int right_uncovered_limit) :
                right_uncovered_limit_(right_uncovered_limit) { }

        VJFilteringInfo Filter(const VJHits &vj_hits) const {
            bool bad = vj_hits.GetJHitByIndex(0).RightUncovered() > right_uncovered_limit_;
            TRACE("Right uncovered (" << vj_hits.GetJHitByIndex(0).RightUncovered() << ") > limit (" <<
                          right_uncovered_limit_ << "): " << bad);
            VJFilteringInfo res(vj_hits.Read());
            if(bad) {
                res.UpdateFilteringReason(VJFilteringReason::RightUncoveredLimitReason);
                res.right_uncovered_length = vj_hits.GetJHitByIndex(0).RightUncovered();
            }
            return res;
        }

    private:
        DECL_LOGGER("RightCoverageFilter");
    };

    class VSegmentLengthFilter : public VJHitsFilter {
        size_t min_v_segment_length_;

    public:
        VSegmentLengthFilter(size_t min_v_segment_length) :
                min_v_segment_length_(min_v_segment_length) { }

        VJFilteringInfo Filter(const VJHits &vj_hits) const {
            bool bad = vj_hits.GetVHitByIndex(0).SegmentLength() < min_v_segment_length_;
            TRACE("V segment length (" << vj_hits.GetVHitByIndex(0).SegmentLength() << ") < limit (" <<
                          min_v_segment_length_ << "): " << bad);
            VJFilteringInfo res(vj_hits.Read());
            if(bad) {
                res.UpdateFilteringReason(VJFilteringReason::VSegmentLengthReason);
                res.v_segment_length = vj_hits.GetVHitByIndex(0).SegmentLength();
            }
            return res;
        }

    private:
        DECL_LOGGER("VSegmentLengthFilter");
    };

    class JSegmentLengthFilter : public VJHitsFilter {
        size_t min_j_segment_length_;

    public:
        JSegmentLengthFilter(size_t min_j_segment_length) :
                min_j_segment_length_(min_j_segment_length) { }

        VJFilteringInfo Filter(const VJHits &vj_hits) const {
            bool bad = vj_hits.GetJHitByIndex(0).SegmentLength() < min_j_segment_length_;
            TRACE("J segment length (" << vj_hits.GetJHitByIndex(0).SegmentLength() << ") < limit (" <<
                  min_j_segment_length_ << "): " << bad);
            VJFilteringInfo res(vj_hits.Read());
            if(bad) {
                res.UpdateFilteringReason(VJFilteringReason::JSegmentLengthReason);
                res.j_segment_length = vj_hits.GetJHitByIndex(0).SegmentLength();
            }
            return res;
        }

    private:
        DECL_LOGGER("JSegmentLengthFilter");
    };

    class AlignedSegmentLengthFilter : public VJHitsFilter {
        size_t min_aligned_segment_length_;

    public:
        AlignedSegmentLengthFilter(size_t min_aligned_segment_length) :
                min_aligned_segment_length_(min_aligned_segment_length) { }

        VJFilteringInfo Filter(const VJHits &vj_hits) const {
            bool bad = vj_hits.AlignedSegmentLength() < min_aligned_segment_length_;
            TRACE("Alignment segment length (" << vj_hits.AlignedSegmentLength() << ") < limit ("
                  << min_aligned_segment_length_ << "): " << bad);
            VJFilteringInfo res(vj_hits.Read());
            if(bad) {
                res.UpdateFilteringReason(VJFilteringReason::AlignedSegmentLengthReason);
                res.aligned_segment_length = vj_hits.AlignedSegmentLength();
            }
            return res;
        }

    private:
        DECL_LOGGER("AlignedSegmentLengthFilter");
    };

    class CustomVjHitsFilter {
        std::vector<std::shared_ptr<VJHitsFilter> > vj_filters_;

    public:
        void Add(std::shared_ptr<VJHitsFilter> vj_filter);

        VJFilteringInfo Filter(const VJHits &vj_hits) const;
    };

    class VersatileVjFilter {
        const VJFinderConfig::AlgorithmParams::FilteringParams &filtering_params_;
        CustomVjHitsFilter custom_vj_filter_;

        void InitializeCustomVjFinder();

    public:
        VersatileVjFilter(const VJFinderConfig::AlgorithmParams::FilteringParams &filtering_params) :
                filtering_params_(filtering_params) {
            InitializeCustomVjFinder();
        }

        VJFilteringInfo Filter(const VJHits &vj_hits) const { return custom_vj_filter_.Filter(vj_hits); }
    };
}