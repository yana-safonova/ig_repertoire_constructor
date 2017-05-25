#pragma once

#include <string>

template<typename T>
T compress_string(const T& str) {
    std::string compressed = "";
    char prev = 0;
    for (char cur : str) {
        if (cur != prev) {
            compressed += cur;
        }
        prev = cur;
    }
    return compressed;
}
