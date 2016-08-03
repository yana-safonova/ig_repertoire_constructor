//
// Created by Andrew Bzikadze on 7/26/16.
//

#pragma once
#include "alignment_utils/alignment_positions.hpp"
#include "vj_alignment_structs.hpp"

namespace vdj_labeler {

class AbstractDAlignmentPositionsCalculator {
public:
    AbstractDAlignmentPositionsCalculator() = default;

    AbstractDAlignmentPositionsCalculator(const AbstractDAlignmentPositionsCalculator &)           = delete;
    AbstractDAlignmentPositionsCalculator& operator=(const AbstractDAlignmentPositionsCalculator&) = delete;
    AbstractDAlignmentPositionsCalculator(AbstractDAlignmentPositionsCalculator &&)                = delete;
    AbstractDAlignmentPositionsCalculator& operator=(AbstractDAlignmentPositionsCalculator&&)      = delete;

    virtual alignment_utils::AlignmentPositions ComputeDAlignmentPositions(
        const std::vector<vj_finder::VGeneHit> &,
        const std::vector<vj_finder::JGeneHit> &) const = 0;

    virtual ~AbstractDAlignmentPositionsCalculator() { }
};

} // End namespace vdj_labeler