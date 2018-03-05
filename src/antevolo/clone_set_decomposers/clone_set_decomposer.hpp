#pragma once

#include <decomposition.hpp>
#include <string>
#include "../clone_set_with_fakes.hpp"

namespace antevolo {
    class CloneSetDecomposer {
    protected:
        const annotation_utils::CDRAnnotatedCloneSet& clone_set_;

        std::string GetGeneBaseName(seqan::CharString name) const;
        virtual std::string GetGeneKeyByClone(const annotation_utils::AnnotatedClone& clone) const = 0;

    public:
        CloneSetDecomposer(const annotation_utils::CDRAnnotatedCloneSet& clone_set) :
                clone_set_(clone_set) { }

        virtual core::Decomposition CreateDecomposition() const = 0;



        virtual  ~CloneSetDecomposer() { }
    };

    typedef std::shared_ptr<CloneSetDecomposer> CloneSetDecomposerPtr;
}