#include "hcdr3_j_labeler.hpp"

namespace cdr_labeler {
    using namespace annotation_utils;

    size_t HCDR3JLabeler::ComputeEndPosition(const germline_utils::ImmuneGene &) {
        VERIFY_MSG(false, "Implement me!");
        return size_t(-1);
    }

    CDRRange HCDR3JLabeler::ComputeRange(const germline_utils::ImmuneGene &immune_gene, CDRRange) {
        size_t end_position = ComputeEndPosition(immune_gene);
        return CDRRange(size_t(-1), end_position);
    }
}