//
// Created by Andrew Bzikadze on 7/26/16.
//

#pragma once

#include "abstract_d_alignment_positions_calculator.hpp"

namespace vdj_labeler {

class CustomDAlignmentPositionsCalculator : public AbstractDAlignmentPositionsCalculator {

    virtual alignment_utils::AlignmentPositions ComputeDAlignmentPositions(
        const std::vector<vj_finder::VGeneHit> &,
        const std::vector<vj_finder::JGeneHit> &) const override;
};

} // End namespace vdj_labeler
