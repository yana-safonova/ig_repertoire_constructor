//
// Created by Andrew Bzikadze on 7/16/16.
//

#pragma once
#include "vdj_config.hpp"
#include "alignment_utils/alignment_positions.hpp"

namespace vdj_labeler {

class AbstractDAlignmentPositionChecker {
public:
    AbstractDAlignmentPositionChecker() = delete;
    AbstractDAlignmentPositionChecker(const AbstractDAlignmentPositionChecker &) = delete;
    AbstractDAlignmentPositionChecker& operator=(const AbstractDAlignmentPositionChecker&) = delete;
    AbstractDAlignmentPositionChecker(AbstractDAlignmentPositionChecker &&) = delete;
    AbstractDAlignmentPositionChecker& operator=(AbstractDAlignmentPositionChecker&&) = delete;

    AbstractDAlignmentPositionChecker(const VDJLabelerConfig::DAlignmentQualityParams) { };
    virtual bool DAlignmentPositionsAreGood(const alignment_utils::AlignmentPositions &d_alignment_positions) const = 0;
};

} // End namespace vdj_labeler
