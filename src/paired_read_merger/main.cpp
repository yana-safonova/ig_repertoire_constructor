#include "reads_merger.hpp"
#include <utils/string_tools.hpp>

#include "logger/log_writers.hpp"

merger_setting parse_settings(int argc, char *argv[]) {
    merger_setting setting;
    string min_overlap_str = "--min-overlap=";
    string max_mismatch_str = "--max-mismatch=";
    string simulated_mode_str = "--simulated-mode";
    for(size_t i = 4; i < static_cast<size_t>(argc); i++) {
        string tmp(argv[i]);
        if(tmp.substr(0, min_overlap_str.size()) == min_overlap_str) {
            tmp = tmp.substr(min_overlap_str.size(), tmp.size() - min_overlap_str.size());
            setting.min_overlap = string_to_number<size_t>(tmp);
        }
        else if(tmp.substr(0, max_mismatch_str.size()) == max_mismatch_str) {
            tmp = tmp.substr(max_mismatch_str.size(), tmp.size() - max_mismatch_str.size());
            setting.max_mismatch_rate = string_to_number<double>(tmp);
        }
        else if(tmp == simulated_mode_str)
            setting.simulated_mode = true;

    }
    setting.print();
    return setting;
}

void create_console_logger(string log_filename) {
    using namespace logging;
    logger *lg = create_logger(path::FileExists(log_filename) ? log_filename : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}


int main(int argc, char *argv[]) {
    /*
     * argv[1] - left fastq reads
     * argv[2] - right fastq reads
     * argv[3] - prefix of output files (prefix.fastq prefix.stats)
     */
    create_console_logger("");

    if(argc < 4) {
        ERROR("paired_read_merger left_reads.fq right_reads.fq output_prefix [--min-overlap=N1 --max-mismatch=N2 --simulated-mode]");
        return 1;
    }

    vector<PairedFastqRead> paired_reads = PairedFastqReader(argv[1],
            argv[2]).Read();
    INFO(paired_reads.size() << " paired reads were read from " << argv[1] <<
            " and " << argv[2]);
    vector<FastqRead> merged_reads = PairedReadsMerger(parse_settings(argc, argv)).Merge(paired_reads);
    INFO(merged_reads.size() << " read from " << paired_reads.size() << " were successfully merged");
    FastqWriter(string(argv[3])).Write(merged_reads);
    INFO("Merged reads were written to " << string(argv[3]));
    return 0;
}
