#include "fastq_read_archive.hpp"

FastqReadArchive::FastqReadArchive(std::string fastq_file_fname) {
    std::vector<CharString> read_headers;
    std::vector<Dna5String> read_seqs;
    seqan::SeqFileIn seqFileIn_reads(fastq_file_fname.c_str());
    readRecords(read_headers, read_seqs, seqFileIn_reads);
    for(size_t i = 0; i < read_seqs.size(); i++)
        reads_.push_back(Read(read_headers[i], read_seqs[i]));
}

size_t FastqReadArchive::size() const {
    return reads_.size();
}


const Read& FastqReadArchive::operator[](size_t index) const {
    assert(index < reads_.size());
    return reads_[index];
}