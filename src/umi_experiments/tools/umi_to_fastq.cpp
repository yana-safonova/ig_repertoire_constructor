#include <ostream>
#include <iostream>
#include <fstream>
#include <string.h>
#include <logger/log_writers.hpp>
#include <seqan/seq_io.h>
#include "umi_utils.hpp"
#include "utils.hpp"

int main(int argc, char* argv[]){
    create_console_logger();
    if((argc != 5 && argc != 7) || strncmp(argv[1], "-i", 3) != 0 || strncmp(argv[3], "-o", 3) != 0 || (argc == 7 && strncmp(argv[5], "-c", 3) != 0)) {
        INFO("Extracts UMIs from fastq file into a separate one.");
        INFO("Usage: -i <input file> -o <output file> [-c <UMI cleavage length>]");
        return 1;
    }
    size_t umi_cleavage_length = (argc == 7) ? std::stoull(std::string(argv[6])) : 0;

    INFO("Reading fastq");
    seqan::SeqFileIn infile(argv[2]);
    std::vector<seqan::CharString> input_ids;
    {
        std::vector<seqan::Dna5String> input_reads;
        std::vector<seqan::DnaQString> input_qual;
        readRecords(input_ids, input_reads, input_qual, infile);
    }
    INFO(input_ids.size() << " records read.");

    INFO("Extracting barcodes.");
    std::vector<seqan::Dna5String> umis;
    std::vector<seqan::DnaQString> umi_quals;
    extract_barcodes_from_read_ids(input_ids, umis, umi_quals);

    INFO("Saving barcodes to " << argv[4]);
    seqan::SeqFileOut outfile(argv[4]);
    for (size_t i = 0; i < umis.size(); i ++) {
        if (umi_quals.empty()) {
            writeRecord(outfile, input_ids[i], seqan::suffix(umis[i], umi_cleavage_length));
        } else {
            writeRecord(outfile, input_ids[i], seqan::suffix(umis[i], umi_cleavage_length), seqan::suffix(umi_quals[i], umi_cleavage_length));
        }
    }
    close(outfile);
    INFO(input_ids.size() << " barcodes were extracted from " << argv[2] << " to " << argv[4]);
}
