#pragma once

#include "aa_annotation.hpp"
#include "../cdr_labeling_primitives.hpp"

namespace annotation_utils {

    class BaseAACalculator {
    public:
        virtual AminoAcidAnnotation<core::Read> ComputeAminoAcidAnnotation(const core::Read& read,
                                                                           const CDRLabeling &cdr_labeling) const = 0;

        virtual ~BaseAACalculator() { }
    };

    class SimpleAACalculator : public BaseAACalculator {
    private:
        bool ComputeInFrame(const CDRLabeling &cdr_labeling) const;

        bool FindStopCodon(const AAString &aa_str) const;

    public:
        AminoAcidAnnotation<core::Read> ComputeAminoAcidAnnotation(const core::Read& read,
                                                                   const CDRLabeling &cdr_labeling) const override;
    };

    using BaseAACalculatorPtr = std::unique_ptr<BaseAACalculator>;
}