#include <assert.h>
#include <ostream>
#include <iostream>
#include <fstream>
#include <string.h>
#include <logger/log_writers.hpp>

#include <seqan/seq_io.h>
#include "../ig_tools/utils/string_tools.hpp"

using seqan::Dna5String;
using seqan::DnaQString;
using seqan::SeqFileIn;
using seqan::SeqFileOut;
using seqan::CharString;

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int argc, char * argv[]){
    create_console_logger();
    if(argc != 5 || strncmp(argv[1], "-i", 3) != 0 || strncmp(argv[3], "-o", 3) != 0) {
        INFO("Extracts UMIs from fastq file into a separate one.");
        INFO("Usage: -i <input file> -o <output file>");
        return 1;
    }

    INFO("Reading fastq");
    SeqFileIn infile(argv[2]);
    std::vector<CharString> input_ids;
    {
        std::vector<Dna5String> input_reads;
        std::vector<DnaQString> input_qual;
        readRecords(input_ids, input_reads, input_qual, infile);
    }

    INFO("Extracting barcodes.");
    SeqFileOut outfile(argv[4]);
    char delim = ' ';
    for (auto& id : input_ids) {
        std::string s = seqan_string_to_string(id);
        std::size_t space = s.find(delim);
        if (space == std::string::npos) {
            std::string delims = " _";
            delim = delims[delims.find(delim) ^ 1];
            space = s.find(delim);
        }
        assert(space != std::string::npos);
        space = s.find(delim, space + 1);
        assert(space != std::string::npos);
        std::string meta = s.substr(0, space);
        std::string umi_info = s.substr(space + 1);
        assert(!umi_info.empty());
        std::size_t colon = umi_info.find(':');
        assert(colon != std::string::npos);
        umi_info = umi_info.substr(colon + 1);
        colon = umi_info.find(':');
        assert(colon != std::string::npos);
        assert(colon * 2 + 1 == umi_info.length());
        auto umi = umi_info.substr(0, colon);
        auto qual = umi_info.substr(colon + 1);
        writeRecord(outfile, meta, umi, qual);
    }

    INFO(input_ids.size() << " barcodes were extracted from " << argv[2] << " to " << argv[4]);
}
