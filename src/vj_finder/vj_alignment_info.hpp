#pragma once

#include "vj_finder_config.hpp"
#include "vj_alignment_structs.hpp"

namespace vj_finder {
    class VJAlignmentInfo {
        std::vector<VJHits> alignment_records_;
        std::vector<core::Read*> filtered_reads_;

    public:
        VJAlignmentInfo() { }

        void UpdateHits(VJHits vj_hits);

        void UpdateFilteredReads(const core::Read& read);

        void Update(VJAlignmentInfo vj_alignment);
    };

    class VJAlignmentOutput {
        const VJFinderConfig::IOParams::OutputParams &output_params_;
        const VJAlignmentInfo &alignment_info_;

        VJAlignmentOutput(const VJFinderConfig::IOParams::OutputParams &output_params,
                          const VJAlignmentInfo &alignment_info) :
                output_params_(output_params), alignment_info_(alignment_info) { }

        void OutputFilteredReads() const;

        void OutputAlignmentInfo() const;

        void OutputCleanedReads() const;
    };
}