//
// Created by Andrew Bzikadze on 3/31/17.
//

#pragma once

#include "verify.hpp"

namespace ig_simulator {

template<class Pointer>
const Pointer& check_pointer(const Pointer& p) {
    VERIFY(p != nullptr);
    return p;
}

} // End namespace ig_simulator
