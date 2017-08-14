#pragma once

#include "clone_set_decomposer.hpp"
#include <string>

namespace antevolo {
    class VJCloneSetDecomposer : public CloneSetDecomposer {
        std::string GetGeneBaseName(seqan::CharString name) const;

        std::string GetVJKeyByClone(const annotation_utils::AnnotatedClone &clone) const;

    public:
        VJCloneSetDecomposer(const annotation_utils::CDRAnnotatedCloneSet& clone_set) :
                CloneSetDecomposer(clone_set) { }

        core::Decomposition CreateDecomposition() const;
    };
}