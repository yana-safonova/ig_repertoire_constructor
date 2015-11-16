#pragma once

#include <seqan/seq_io.h>

using seqan::Dna5String;
using seqan::CharString;
using seqan::length;
using seqan::SeqFileIn;
using seqan::SeqFileOut;

struct Read {
    CharString name;
    Dna5String seq;

    Read(CharString new_name, Dna5String new_seq) :
            name(new_name), seq(new_seq) { }
};

class FastqReadArchive {
    std::vector<Read> reads_;

public:
    FastqReadArchive(std::string fastq_file_fname);

    size_t size() const;

    typedef std::vector<Read>::const_iterator read_iterator;

    read_iterator cbegin() const { return reads_.cbegin(); }

    read_iterator cend() const { return reads_.cend(); }

    const Read& operator[](size_t index) const;

};