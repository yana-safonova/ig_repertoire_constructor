#pragma once

#include "alignment_utils/pairwise_alignment.hpp"

namespace vdj_labeler {

class AlignmentQualityChecker {
public:
    AlignmentQualityChecker() { }
    AlignmentQualityChecker(const AlignmentQualityChecker &) = delete;
    AlignmentQualityChecker& operator=(const AlignmentQualityChecker&) = delete;
    AlignmentQualityChecker(AlignmentQualityChecker &&) = delete;
    AlignmentQualityChecker& operator=(AlignmentQualityChecker&&) = delete;

    virtual bool AlignmentIsGood(const alignment_utils::ImmuneGeneReadAlignment &ig_gene_alignment) const = 0;
    virtual ~AlignmentQualityChecker() { }
};

} // End namespace vdj_labeler