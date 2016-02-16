#include "pairing_fastq_utils.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string.hpp>

#include <sstream>

const std::string HeaderConfig::isotype_prefix = "PRCONS=";

const std::string HeaderConfig::main_delimeters = "|";

const std::string HeaderConfig::barcode_delimeters = "_";

const std::string HeaderConfig::presize_symbol = "=";

const std::string HeaderConfig::size_delimeter = ":";

bool PairingFastqUtils::LineIsHeader(std::string line) {
    if(line[0] != '@')
        return false;
    std::vector<std::string> splits;
    boost::split(splits, line, boost::is_any_of(HeaderConfig::main_delimeters));
    return splits.size() == HeaderConfig::num_header_splits;
}

bool PairingFastqUtils::LineIsPlus(std::string line) {
    return line == "+";
}

std::string PairingFastqUtils::ExtractDropletBarcode(std::string header) {
    assert(header[0] == '@');
    size_t index = header.find(HeaderConfig::barcode_delimeters);
    return header.substr(1, index - 1);
}

IgIsotype PairingFastqUtils::ExtractIsotypeFromHeader(std::string header) {
    std::vector<std::string> splits;
    boost::split(splits, header, boost::is_any_of(HeaderConfig::main_delimeters));
    assert(splits.size() == HeaderConfig::num_header_splits);
    std::string isotype_str = splits[HeaderConfig::isotype_index].substr(HeaderConfig::isotype_prefix.size());
    return IgIsotype(isotype_str);
}

std::string PairingFastqUtils::ExtractUmiFromHeader(std::string header) {
    size_t start_index = header.find(HeaderConfig::barcode_delimeters) + 1;
    size_t end_index = header.find(HeaderConfig::main_delimeters);
    return header.substr(start_index, end_index - start_index);
}

size_t PairingFastqUtils::ExtractSizeFromHeader(std::string header) {
    size_t start_index = header.find(HeaderConfig::presize_symbol) + 1;
    size_t end_index = header.find(HeaderConfig::size_delimeter);
    size_t mb_size;
    std::stringstream ss;
    ss << header.substr(start_index, end_index - start_index);
    ss >> mb_size;
    return mb_size;
}

seqan::Dna5String PairingFastqUtils::ConvertToDnaString(std::string sequence) {
    std::replace(sequence.begin(), sequence.end(), '-', 'N');
    return seqan::Dna5String(sequence);
}