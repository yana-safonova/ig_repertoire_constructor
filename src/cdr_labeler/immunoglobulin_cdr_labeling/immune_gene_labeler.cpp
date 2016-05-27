#include "immune_gene_labeler.hpp"

namespace cdr_labeler {
    CDRLabeling ImmuneGeneCDRLabeler::ComputeLabeling(const germline_utils::ImmuneGene &immune_gene) {
        VERIFY_MSG(immune_gene.GeneType() == gene_type_, "Type of immune gene (" << immune_gene.GeneType() <<
                                                         ") is not " << gene_type_);
        CDRRange cdr1 = cdr1_labeler_ptr_->ComputeLoopRange(immune_gene, CDRRange());
        CDRRange cdr2 = cdr2_labeler_ptr_->ComputeLoopRange(immune_gene, cdr1);
        CDRRange cdr3 = cdr3_labeler_ptr_->ComputeLoopRange(immune_gene, cdr2);
        return CDRLabeling(cdr1, cdr2, cdr3);
    }

    SingleLoopLabelerPtr ImmuneGeneCDRHelper::GetCDR1Labeler(germline_utils::ImmuneGeneType gene_type,
                                                             const CDRLabelerConfig::CDRsParams &cdr_params) {
        if(gene_type.Segment() == germline_utils::SegmentType::JoinSegment)
            return SingleLoopLabelerPtr(new SingleLoopLabeler(gene_type));
        return SingleLoopLabelerPtr(new HCDR1Labeler(cdr_params.hcdr1_params));
    }

    SingleLoopLabelerPtr ImmuneGeneCDRHelper::GetCDR2Labeler(germline_utils::ImmuneGeneType gene_type,
                                                             const CDRLabelerConfig::CDRsParams &cdr_params) {
        if(gene_type.Segment() == germline_utils::SegmentType::JoinSegment)
            return SingleLoopLabelerPtr(new SingleLoopLabeler(gene_type));
        return SingleLoopLabelerPtr(new HCDR2Labeler(cdr_params.hcdr2_params));
    }

    SingleLoopLabelerPtr ImmuneGeneCDRHelper::GetCDR3Labeler(germline_utils::ImmuneGeneType gene_type,
                                                             const CDRLabelerConfig::CDRsParams &cdr_params) {
        if(gene_type.Segment() == germline_utils::SegmentType::JoinSegment)
            return SingleLoopLabelerPtr(new HCDR3JLabeler(cdr_params.hcdr3_params));
        return SingleLoopLabelerPtr(new HCDR3VLabeler(cdr_params.hcdr3_params));
    }
}