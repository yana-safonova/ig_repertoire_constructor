#include "hcdr3_v_labeler.hpp"

namespace cdr_labeler {
    using namespace annotation_utils;

    size_t HCDR3VLabeler::ComputeStartPosition(const germline_utils::ImmuneGene &, CDRRange) {
        return size_t(-1);
    }

    CDRRange HCDR3VLabeler::ComputeRange(const germline_utils::ImmuneGene &immune_gene, CDRRange previous_cdr) {
        //size_t start_pos = ComputeStartPosition(immune_gene, previous_cdr);
        return CDRRange(ComputeStartPosition(immune_gene, previous_cdr), size_t(-1));
    }
}