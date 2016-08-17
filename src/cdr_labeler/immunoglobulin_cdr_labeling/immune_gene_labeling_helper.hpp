#pragma once

#include "immune_gene_labeler.hpp"
#include "annotated_gene_labeler.hpp"

namespace cdr_labeler {
    class ImmuneGeneLabelingHelper {
    public:
        static BaseImmuneGeneCDRLabelerPtr GetAnnotatedLabeler(germline_utils::SegmentType segment_type,
                                                               const CDRLabelerConfig::CDRsParams &cdr_params);
    };
}