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

template <class T>
T check_numeric_nonnegative(const T x) {
    VERIFY(x >= 0);
    return x;
}

template <class T>
T check_numeric_positive(const T x) {
    VERIFY(x > 0);
    return x;
}

} // End namespace ig_simulator
