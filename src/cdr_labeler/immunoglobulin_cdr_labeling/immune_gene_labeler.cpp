#include <verify.hpp>

#include "immune_gene_labeler.hpp"

namespace cdr_labeler {
    using namespace annotation_utils;

    CDRLabeling DeNovoImmuneGeneCDRLabeler::ComputeLabeling(const germline_utils::ImmuneGene &immune_gene) {
        VERIFY_MSG(immune_gene.GeneType() == gene_type_, "Type of immune gene (" << immune_gene.GeneType() <<
                                                         ") is not " << gene_type_);
        CDRRange cdr1 = cdr1_labeler_ptr_->ComputeLoopRange(immune_gene, CDRRange());
        CDRRange cdr2 = cdr2_labeler_ptr_->ComputeLoopRange(immune_gene, cdr1);
        CDRRange cdr3 = cdr3_labeler_ptr_->ComputeLoopRange(immune_gene, cdr2);
        return CDRLabeling(cdr1, cdr2, cdr3);
    }
}