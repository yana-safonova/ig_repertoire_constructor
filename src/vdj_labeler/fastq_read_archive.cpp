#include "fastq_read_archive.hpp"

#include "seqan/sequence.h"

using namespace std;

std::ostream& operator<<(std::ostream& out, const Read& read) {
    out << "Name: " << read.name << ", seq: " << read.seq;
    return out;
}

//-----------------------------------------------------------------

FastqReadArchive::FastqReadArchive(std::string fastq_file_fname) {
    std::vector<CharString> read_headers;
    std::vector<Dna5String> read_seqs;
    seqan::SeqFileIn seqFileIn_reads(fastq_file_fname.c_str());
    readRecords(read_headers, read_seqs, seqFileIn_reads);
    for(size_t i = 0; i < read_seqs.size(); i++) {
        ReadPtr read_ptr(new Read(read_headers[i], read_seqs[i]));
        reads_.push_back(read_ptr);
        name_read_map_[string(toCString(read_headers[i]))] = read_ptr;
    }
}

size_t FastqReadArchive::size() const {
    return reads_.size();
}


ReadPtr FastqReadArchive::operator[](size_t index) const {
    assert(index < reads_.size());
    return reads_[index];
}


ReadPtr FastqReadArchive::GetReadByName(std::string read_name) const {
    assert(name_read_map_.find(read_name) != name_read_map_.end());
    return name_read_map_.at(read_name);
}