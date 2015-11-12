#include "fastq_read_archive.hpp"

FastqReadArchive::FastqReadArchive(std::string fastq_file_fname) {
    seqan::SeqFileIn seqFileIn_reads(fastq_file_fname.c_str());
    readRecords(read_headers_, read_seqs_, seqFileIn_reads);
}

size_t FastqReadArchive::size() {
    return read_headers_.size();
}