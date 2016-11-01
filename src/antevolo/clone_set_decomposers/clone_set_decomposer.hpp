#pragma once

#include <decomposition.hpp>
#include "../clone_set_with_fakes.hpp"

namespace antevolo {
    class CloneSetDecomposer {
    protected:
        const CloneSetWithFakes& clone_set_;

    public:
        CloneSetDecomposer(const CloneSetWithFakes& clone_set) :
                clone_set_(clone_set) { }

        virtual core::Decomposition CreateDecomposition() const = 0;

        virtual  ~CloneSetDecomposer() { }
    };
}