//
// Created by Andrew Bzikadze on 7/17/16.
//

#include "abstract_d_alignment_position_checker.hpp"

namespace vdj_labeler {
class InfoBasedDAlignmentPositionChecker : public AbstractDAlignmentPositionChecker {
protected:
    unsigned min_coverage_;

public:
    InfoBasedDAlignmentPositionChecker(const VDJLabelerConfig::DAlignmentQualityParams config_) :
        AbstractDAlignmentPositionChecker(config_),
        min_coverage_(config_.min_coverage)
    { }

    bool DAlignmentPositionsAreGood(const alignment_utils::AlignmentPositions &d_alignment_positions) const override {
        return d_alignment_positions.QueryAlignmentLength() >= min_coverage_;
    }
};
} // End namespace vdj_labeler
