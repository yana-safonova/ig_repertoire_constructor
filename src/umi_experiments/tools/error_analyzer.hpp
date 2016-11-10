#pragma once

#include <string>
#include <seqan/seq_io.h>

class ErrorAnalyzer {
public:
    void readData(std::string input_file_path);
    void performAnalysis();

private:
    std::vector<seqan::CharString> ids_;
    std::vector<seqan::Dna5String> reads_;
};
