//
// Created by Andrew Bzikadze on 3/27/17.
//

#pragma once

#include "metaroot/metaroot.hpp"
#include "annotation_utils/aa_annotation/aa_calculator.hpp"

namespace ig_simulator {

class ProductivityChecker {
private:
    const annotation_utils::BaseAACalculatorPtr aa_calculator;

public:
    ProductivityChecker(annotation_utils::BaseAACalculatorPtr aa_calculator):
        aa_calculator(std::move(aa_calculator))
    { }

    bool IsProductive(const AbstractMetaroot& root) {
        if (root.CDRLabeling().Empty())
            return false;
        core::Read read("", root.Sequence(), 0);
        auto aa = aa_calculator->ComputeAminoAcidAnnotation(read, root.CDRLabeling());
        return not aa.HasStopCodon() and aa.InFrame();
    }

    bool IsProductive(const AbstractMetarootCPtr& root) {
        return IsProductive(*check_pointer(root));
    }
};

} // End namespace ig_simulator
