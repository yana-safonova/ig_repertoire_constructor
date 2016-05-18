#pragma once

#include "vj_finder_config.hpp"
#include "vj_alignment_structs.hpp"

namespace vj_finder {
    class VJAlignmentInfo {
        std::vector<VJHits> alignment_records_;
        std::vector<const core::Read*> filtered_reads_;

    public:
        VJAlignmentInfo() { }

        void UpdateHits(VJHits vj_hits);

        void UpdateFilteredReads(const core::Read& read);

        void Update(VJAlignmentInfo vj_alignment_info);

        size_t NumFilteredReads() const { return filtered_reads_.size(); }

        size_t NumVJHits() const { return alignment_records_.size(); }

        const VJHits& GetVJHitsByIndex(size_t index) const;

        const core::Read& GetFilteredReadByIndex(size_t index) const;
    };

    class VJAlignmentOutput {
        const VJFinderConfig::IOParams::OutputParams &output_params_;
        const VJAlignmentInfo &alignment_info_;

    public:
        VJAlignmentOutput(const VJFinderConfig::IOParams::OutputParams &output_params,
                          const VJAlignmentInfo &alignment_info) :
                output_params_(output_params), alignment_info_(alignment_info) { }

        void OutputFilteredReads() const;

        void OutputAlignmentInfo() const;

        void OutputCleanedReads() const;
    };
}