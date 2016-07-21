#pragma once

#include <seqan/sequence.h>

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

    std::string dna5String_to_string(seqan::Dna5String);

    template<typename SeqanStringType>
    std::string seqan_string_to_string(SeqanStringType seqan_str) {
        std::stringstream ss;
        ss << seqan_str;
        return ss.str();
    }
}