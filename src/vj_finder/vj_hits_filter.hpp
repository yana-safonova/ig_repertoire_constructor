#pragma once

#include "vj_finder_config.hpp"
#include "vj_alignment_structs.hpp"

namespace vj_finder {
    class VJHitsFilter {
    public:
        // returns true if VJ Hits are bad, otherwise returns false
        virtual bool Filter(const VJHits &vj_hits) const = 0;

        ~VJHitsFilter() { }
    };

    class EmptyHitsFilter : public VJHitsFilter {
    public:
        bool Filter(const VJHits &vj_hits) const {
            bool empty = vj_hits.NumJHits() == 0 or vj_hits.NumVHits() == 0;
            TRACE("V hits are empty and J hits are empty: " << empty);
            return vj_hits.NumJHits() == 0 or vj_hits.NumVHits() == 0;
        }

    private:
        DECL_LOGGER("EmptyHitsFilter");
    };

    class LeftCoverageFilter : public VJHitsFilter {
        int left_uncovered_limit_;

    public:
        LeftCoverageFilter(int left_uncovered_limit) :
                left_uncovered_limit_(left_uncovered_limit) { }

        bool Filter(const VJHits &vj_hits) const {
            bool bad = vj_hits.GetVHitByIndex(0).LeftUncovered() > left_uncovered_limit_;
            TRACE("Left uncovered (" << vj_hits.GetVHitByIndex(0).LeftUncovered() <<
                    ") > limit (" << left_uncovered_limit_ << "): " << bad);
            return vj_hits.GetVHitByIndex(0).LeftUncovered() > left_uncovered_limit_;
        }

    private:
        DECL_LOGGER("LeftCoverageFilter");
    };

    class RightCoverageFilter : public VJHitsFilter {
        int right_uncovered_limit_;

    public:
        RightCoverageFilter(int right_uncovered_limit) :
                right_uncovered_limit_(right_uncovered_limit) { }

        bool Filter(const VJHits &vj_hits) const {
            bool bad = vj_hits.GetJHitByIndex(0).RightUncovered() > right_uncovered_limit_;
            TRACE("Right uncovered (" << vj_hits.GetJHitByIndex(0).RightUncovered() << ") > limit (" <<
                          right_uncovered_limit_ << "): " << bad);
            return vj_hits.GetJHitByIndex(0).RightUncovered() > right_uncovered_limit_;
        }

    private:
        DECL_LOGGER("RightCoverageFilter");
    };

    class VSegmentLengthFilter : public VJHitsFilter {
        size_t min_v_segment_length_;

    public:
        VSegmentLengthFilter(size_t min_v_segment_length) :
                min_v_segment_length_(min_v_segment_length) { }

        bool Filter(const VJHits &vj_hits) const {
            bool bad = vj_hits.GetVHitByIndex(0).SegmentLength() < min_v_segment_length_;
            TRACE("V segment length (" << vj_hits.GetVHitByIndex(0).SegmentLength() << ") < limit (" <<
                          min_v_segment_length_ << "): " << bad);
            return vj_hits.GetVHitByIndex(0).SegmentLength() < min_v_segment_length_;
        }

    private:
        DECL_LOGGER("VSegmentLengthFilter");
    };

    class JSegmentLengthFilter : public VJHitsFilter {
        size_t min_j_segment_length_;

    public:
        JSegmentLengthFilter(size_t min_j_segment_length) :
                min_j_segment_length_(min_j_segment_length) { }

        bool Filter(const VJHits &vj_hits) const {
            bool bad = vj_hits.GetJHitByIndex(0).SegmentLength() < min_j_segment_length_;
            TRACE("J segment length (" << vj_hits.GetJHitByIndex(0).SegmentLength() << ") < limit (" <<
                  min_j_segment_length_ << "): " << bad);
            return vj_hits.GetJHitByIndex(0).SegmentLength() < min_j_segment_length_;
        }

    private:
        DECL_LOGGER("JSegmentLengthFilter");
    };

    class AlignedSegmentLengthFilter : public VJHitsFilter {
        size_t min_aligned_segment_length_;

    public:
        AlignedSegmentLengthFilter(size_t min_aligned_segment_length) :
                min_aligned_segment_length_(min_aligned_segment_length) { }

        bool Filter(const VJHits &vj_hits) const {
            bool bad = vj_hits.AlignedSegmentLength() < min_aligned_segment_length_;
            TRACE("Alignment segment length (" << vj_hits.AlignedSegmentLength() << ") < limit ("
                  << min_aligned_segment_length_ << "): " << bad);
            return vj_hits.AlignedSegmentLength() < min_aligned_segment_length_;
        }

    private:
        DECL_LOGGER("AlignedSegmentLengthFilter");
    };

    class CustomVjHitsFilter {
        std::vector<std::shared_ptr<VJHitsFilter> > vj_filters_;

    public:
        void Add(std::shared_ptr<VJHitsFilter> vj_filter);

        bool Filter(const VJHits &vj_hits) const;
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

        bool Filter(const VJHits &vj_hits) const { return custom_vj_filter_.Filter(vj_hits); }
    };
}