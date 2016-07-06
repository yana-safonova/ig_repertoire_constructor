#pragma once

#include "alignment_utils/pairwise_alignment.hpp"

namespace vdj_labeler {

class AlignmentQualityController {
public:
    virtual bool AlignmentIsGood(alignment_utils::ImmuneGeneReadAlignmentPtr ig_gene_alignment) = 0;
    virtual ~AlignmentEstimator() { }
};

} // End namespace vdj_labeler