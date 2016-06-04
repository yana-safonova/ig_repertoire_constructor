#pragma once

#include "cdr_config.hpp"
#include <annotation_utils/annotated_clone_set.hpp>
#include <vj_alignment_info.hpp>

namespace cdr_labeler {
    class CDRLabelingWriter {
        const CDRLabelerConfig::OutputParams &output_config_;
        const vj_finder::VJAlignmentInfo &alignment_info_;
        const annotation_utils::CDRAnnotatedCloneSet &clone_set_;

    public:
        CDRLabelingWriter(const CDRLabelerConfig::OutputParams &output_config,
                          const vj_finder::VJAlignmentInfo &alignment_info,
                          const annotation_utils::CDRAnnotatedCloneSet &clone_set) : output_config_(output_config),
                                                                                     alignment_info_(alignment_info),
                                                                                     clone_set_(clone_set) { }

        void OutputCDRDetails() const;

        void OutputCDR3Fasta() const;
    };
}