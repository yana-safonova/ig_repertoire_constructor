//
// Created by Andrew Bzikadze on 3/27/17.
//

#pragma once

#include "base_repertoire/metaroot/metaroot.hpp"
#include "annotation_utils/aa_annotation/aa_calculator.hpp"

namespace ig_simulator {

class ProductivityChecker {
private:
    const annotation_utils::BaseAACalculatorPtr aa_calculator;

public:
    explicit ProductivityChecker(annotation_utils::BaseAACalculatorPtr aa_calculator =
                        annotation_utils::BaseAACalculatorPtr(new annotation_utils::SimpleAACalculator())):
        aa_calculator(std::move(aa_calculator))
    { }

    bool IsProductive(const AbstractMetaroot& root) const;

    bool IsProductive(const AbstractMetarootCPtr& root) const {
        return IsProductive(*check_pointer(root));
    }

    ProductivityChecker(const ProductivityChecker&) = delete;
    ProductivityChecker(ProductivityChecker&&) = delete;
    ProductivityChecker& operator=(const ProductivityChecker&) = delete;
    ProductivityChecker& operator=(ProductivityChecker&&) = delete;
};

} // End namespace ig_simulator
