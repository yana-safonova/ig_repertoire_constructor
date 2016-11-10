#include <seqan/seq_io.h>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
#include <logger/logger.hpp>

bool create_new_directory(boost::filesystem::path path) {
    namespace fs = boost::filesystem;
    bool exists = fs::exists(path);
    if (exists) {
        fs::remove_all(path);
    }
    fs::create_directory(path);
    return !exists;
}

void write_seqan_records(boost::filesystem::path path, std::vector<seqan::CharString> ids, std::vector<seqan::Dna5String> reads) {
    namespace fs = boost::filesystem;
    seqan::SeqFileOut clusters_file(path.string().c_str());
    writeRecords(clusters_file, ids, reads);
}

void read_seqan_records(const std::string& input_file_name, std::vector<seqan::CharString>& ids, std::vector<seqan::Dna5String>& reads) {
    INFO("Reading sequences from " << input_file_name);
    seqan::SeqFileIn reads_file(input_file_name.c_str());
    readRecords(ids, reads, reads_file);
    INFO(reads.size() << " sequences read");
}
