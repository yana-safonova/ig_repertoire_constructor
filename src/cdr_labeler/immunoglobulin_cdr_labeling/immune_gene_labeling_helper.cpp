#include "immune_gene_labeling_helper.hpp"

namespace cdr_labeler {
    BaseImmuneGeneCDRLabelerPtr ImmuneGeneLabelingHelper::GetAnnotatedLabeler(germline_utils::SegmentType segment_type,
                                                                              const CDRLabelerConfig::CDRsParams &cdr_params) {
        if(segment_type == germline_utils::SegmentType::VariableSegment)
            return BaseImmuneGeneCDRLabelerPtr(new AnnotatedVGeneCDRLabeler(cdr_params.annotated_search_params));
        if(segment_type == germline_utils::SegmentType::JoinSegment)
            return BaseImmuneGeneCDRLabelerPtr(new AnnotatedJGeneCDRLabeler(cdr_params.annotated_search_params));
        VERIFY_MSG(false, "Segment type " << segment_type << " is not variable or join");
        return BaseImmuneGeneCDRLabelerPtr(new BaseImmuneGeneCDRLabeler());
    }
}