#include <verify.hpp>

#include "single_loop_labeler.hpp"

namespace cdr_labeler {
    using namespace annotation_utils;

    void SingleLoopLabeler::CheckConsistencyFatal(const germline_utils::ImmuneGene &immune_gene) {
        VERIFY_MSG(immune_gene.GeneType() == gene_type_, "Type of immune gene (" << immune_gene.GeneType() <<
                                                         ") is not " << gene_type_);
    }

    CDRRange SingleLoopLabeler::ComputeLoopRange(const germline_utils::ImmuneGene &immune_gene, CDRRange previous_cdr) {
        CheckConsistencyFatal(immune_gene);
        return ComputeRange(immune_gene, previous_cdr);
    }
}