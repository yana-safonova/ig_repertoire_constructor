#pragma once
#include "clone_set_decomposer.hpp"

namespace antevolo {
    class VCloneSetDecomposer : public CloneSetDecomposer {
    public:
        VCloneSetDecomposer(const annotation_utils::CDRAnnotatedCloneSet& clone_set) :
                CloneSetDecomposer(clone_set) { }

        std::string GetGeneKeyByClone(const annotation_utils::AnnotatedClone& clone) const;
        core::Decomposition CreateDecomposition() const;

    };
}

