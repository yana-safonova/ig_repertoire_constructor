#pragma once

#include "cdr_config.hpp"
#include <annotation_utils/annotated_clone_set.hpp>

namespace cdr_labeler {
    class CDRLabelingWriter {
        const CDRLabelerConfig::OutputParams &output_config_;
        const annotation_utils::CDRAnnotatedCloneSet &clone_set_;

    public:
        CDRLabelingWriter(const CDRLabelerConfig::OutputParams &output_config,
                          const annotation_utils::CDRAnnotatedCloneSet &clone_set) : output_config_(output_config),
                                                                                     clone_set_(clone_set) { }

        void OutputCDRDetails() const;

        void OutputCDR3Fasta() const;
    };
}