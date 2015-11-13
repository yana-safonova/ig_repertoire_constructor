#pragma once

#include <seqan/seq_io.h>

using seqan::Dna5String;
using seqan::CharString;
using seqan::length;
using seqan::SeqFileIn;
using seqan::SeqFileOut;

class FastqReadArchive {
    std::vector<CharString> read_headers_;
    std::vector<Dna5String> read_seqs_;

public:
    FastqReadArchive(std::string fastq_file_fname);

    size_t size();

};