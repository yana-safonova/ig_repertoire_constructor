#include "immune_gene_labeler.hpp"

namespace cdr_labeler {
    CDRLabeling ImmuneGeneCDRLabeler::ComputeLabeling(const germline_utils::ImmuneGene &immune_gene) {
        VERIFY_MSG(immune_gene.GeneType() == gene_type_, "Type of immune gene (" << immune_gene.GeneType() <<
                                                         ") is not " << gene_type_);
        CDRRange cdr1 = cdr1_labeler_ptr_->ComputeLoopRange(immune_gene, CDRRange());
        CDRRange cdr2 = cdr2_labeler_ptr_->ComputeLoopRange(immune_gene, cdr1);
        CDRRange cdr3 = cdr1_labeler_ptr_->ComputeLoopRange(immune_gene, cdr2); // change into 3!
        return CDRLabeling(cdr1, cdr2, cdr3);
    }

    std::shared_ptr<SingleLoopLabeler> ImmuneGeneCDRHelper::GetCDR1Labeler(germline_utils::ImmuneGeneType gene_type,
                                                                           const CDRLabelerConfig::CDRsParams &cdr_params) {
        return std::shared_ptr<SingleLoopLabeler>(new HCDR1Labeler(cdr_params.hcdr1_params));
    }

    std::shared_ptr<SingleLoopLabeler> ImmuneGeneCDRHelper::GetCDR2Labeler(germline_utils::ImmuneGeneType gene_type,
                                                                           const CDRLabelerConfig::CDRsParams &cdr_params) {
        return std::shared_ptr<SingleLoopLabeler>(new HCDR2Labeler(cdr_params.hcdr2_params));
    }
}