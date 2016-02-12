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
    int qual_mask = 0;
    for (auto& id : input_ids) {
        std::string s = seqan_string_to_string(id);
        auto split_by_umi = split(s, "UMI");
        VERIFY_MSG(split_by_umi.size() <= 2, "Too much 'UMI' strings in read id");
        if (split_by_umi.size() == 1) {
            split_by_umi = split(s, "BARCODE");
            VERIFY_MSG(split_by_umi.size() > 1, "Could not find both 'UMI' and 'BARCODE' in read id");
            VERIFY_MSG(split_by_umi.size() < 3, "Too much 'BARCODE' strings in read id");
        }
        std::string meta = trim(split_by_umi[0]);
        std::string umi_info = split_by_umi[1];
        assert(!umi_info.empty());
        std::size_t colon = umi_info.find(':');
        assert(colon != std::string::npos);
        umi_info = umi_info.substr(colon + 1);
        colon = umi_info.find(':');
        if (colon == std::string::npos) {
            qual_mask |= 1;
            writeRecord(outfile, meta, umi_info);
        } else {
            qual_mask |= 2;
            assert(colon * 2 + 1 == umi_info.length());
            auto umi = umi_info.substr(0, colon);
            auto qual = umi_info.substr(colon + 1);
            writeRecord(outfile, meta, umi, qual);
        }
    }
    close(outfile);
    if (qual_mask == 3) {
        ERROR("Found both UMIs with quality data and without.");
    }
    INFO(input_ids.size() << " barcodes were extracted from " << argv[2] << " to " << argv[4]);
}
