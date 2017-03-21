#include "convert.hpp"

#include <seqan/stream.h>

namespace core {
    std::string dna5String_to_string(seqan::Dna5String dna_string) {
        std::stringstream ss;
        ss << dna_string;
        return ss.str();
    }
}