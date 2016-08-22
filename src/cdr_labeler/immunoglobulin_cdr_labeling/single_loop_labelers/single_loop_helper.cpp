#include "single_loop_helper.hpp"

namespace cdr_labeler {
    SingleLoopLabelerPtr SingleLoopLabelingHelper::GetCDR1Labeler(germline_utils::ImmuneGeneType gene_type,
                                                                  const CDRLabelerConfig::CDRsParams &cdr_params) {
        if(gene_type.Segment() == germline_utils::SegmentType::JoinSegment)
            return SingleLoopLabelerPtr(new SingleLoopLabeler(gene_type));
        return SingleLoopLabelerPtr(new HCDR1Labeler(cdr_params.hcdr1_params));
    }

    SingleLoopLabelerPtr SingleLoopLabelingHelper::GetCDR2Labeler(germline_utils::ImmuneGeneType gene_type,
                                                                  const CDRLabelerConfig::CDRsParams &cdr_params) {
        if(gene_type.Segment() == germline_utils::SegmentType::JoinSegment)
            return SingleLoopLabelerPtr(new SingleLoopLabeler(gene_type));
        return SingleLoopLabelerPtr(new HCDR2Labeler(cdr_params.hcdr2_params));
    }

    SingleLoopLabelerPtr SingleLoopLabelingHelper::GetCDR3Labeler(germline_utils::ImmuneGeneType gene_type,
                                                                  const CDRLabelerConfig::CDRsParams &cdr_params) {
        if(gene_type.Segment() == germline_utils::SegmentType::JoinSegment)
            return SingleLoopLabelerPtr(new HCDR3JLabeler(cdr_params.hcdr3_params));
        return SingleLoopLabelerPtr(new HCDR3VLabeler(cdr_params.hcdr3_params));
    }
}