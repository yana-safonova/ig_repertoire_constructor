#pragma once

#include "../clone_set_with_fakes.hpp"
namespace antevolo {
    class BaseCandidateCalculator {
    protected:
        CloneSetWithFakesPtr clone_set_ptr_;

    public:
        BaseCandidateCalculator(CloneSetWithFakesPtr clone_set_ptr) :
                clone_set_ptr_(clone_set_ptr) { }

        //virtual std::vector<ClonallyRelatedCandidates> ComputeCandidates(core::DecompositionClass decomposition_class) = 0;

        virtual ~BaseCandidateCalculator() { }
    };
}