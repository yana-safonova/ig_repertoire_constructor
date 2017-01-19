#pragma once

#include "../clone_set_with_fakes.hpp"
namespace antevolo {
    class BaseCandidateCalculator {
    protected:
        const annotation_utils::CDRAnnotatedCloneSet& clone_set_;

    public:
        BaseCandidateCalculator(const annotation_utils::CDRAnnotatedCloneSet& clone_set) : clone_set_(clone_set) { }

        //virtual std::vector<ClonallyRelatedCandidates> ComputeCandidates(core::DecompositionClass decomposition_class) = 0;

        virtual ~BaseCandidateCalculator() { }
    };
}