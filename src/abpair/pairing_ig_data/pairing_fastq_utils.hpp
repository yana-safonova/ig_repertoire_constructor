#pragma once

#include "ig_isotype.hpp"

#include <string>
#include "seqan/sequence.h"

struct HeaderConfig {
    static const std::string isotype_prefix;
    static const size_t num_header_splits = 3;
    static const size_t isotype_index = 2;
    static const std::string main_delimeters;
    static const std::string barcode_delimeters;
    static const std::string presize_symbol;
    static const std::string size_delimeter;
};

class PairingFastqUtils {
public:
    static bool LineIsHeader(std::string line);

    static bool LineIsPlus(std::string line);

    static std::string ExtractDropletBarcode(std::string header);

    static IgIsotype ExtractIsotypeFromHeader(std::string header);

    static std::string ExtractUmiFromHeader(std::string header);

    static size_t ExtractSizeFromHeader(std::string);

    static seqan::Dna5String ConvertToDnaString(std::string sequence);
};