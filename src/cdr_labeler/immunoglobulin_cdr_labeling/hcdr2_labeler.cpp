#include "hcdr2_labeler.hpp"

namespace cdr_labeler {
    size_t HCDR2Labeler::ComputeStartPosition(const germline_utils::ImmuneGene &immune_gene,
                                              CDRRange previous_cdr) {
        VERIFY_MSG(false, "Implement me!");
        return size_t(-1);
    }

    size_t HCDR2Labeler::ComputeEndPosition(const germline_utils::ImmuneGene &immune_gene, size_t start_pos) {
        VERIFY_MSG(false, "Implement me!");
        return size_t(-1);
    }

    CDRRange HCDR2Labeler::ComputeRange(const germline_utils::ImmuneGene &immune_gene, CDRRange previous_cdr) {
        size_t start_pos = ComputeStartPosition(immune_gene, previous_cdr);
        size_t end_pos = ComputeEndPosition(immune_gene, start_pos);
        return CDRRange(start_pos, end_pos);
    }
}