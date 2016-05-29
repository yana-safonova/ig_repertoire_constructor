#pragma once

#include <sstream>

namespace core {
    template<typename T1, typename T2>
    T2 convert(T1 obj) {
        std::stringstream ss;
        ss << obj;
        T2 converted;
        ss >> converted;
        return converted;
    };
}