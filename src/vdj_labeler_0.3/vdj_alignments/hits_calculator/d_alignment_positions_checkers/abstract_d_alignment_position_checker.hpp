//
// Created by Andrew Bzikadze on 7/16/16.
//

#pragma once
#include "vdj_config.hpp"
#include "alignment_utils/alignment_positions.hpp"

namespace vdj_labeler {

class AbstractDAlignmentPositionChecker {
public:
    AbstractDAlignmentPositionChecker(const VDJLabelerConfig::DAlignmentQualityParams) { };
    virtual bool DAlignmentPositionsAreGood(const alignment_utils::AlignmentPositions &d_alignment_positions) const = 0;
};

} // End namespace vdj_labeler
