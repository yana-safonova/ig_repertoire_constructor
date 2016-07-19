#pragma once

#include "alignment_utils/pairwise_alignment.hpp"

namespace vdj_labeler {

class AlignmentQualityChecker {
public:
    virtual bool AlignmentIsGood(alignment_utils::ImmuneGeneReadAlignmentPtr ig_gene_alignment) const = 0;
    virtual ~AlignmentQualityChecker() { }
};

} // End namespace vdj_labeler