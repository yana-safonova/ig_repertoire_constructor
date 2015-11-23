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

    Read() : name(), seq() { }
};

std::ostream& operator<<(std::ostream& out, const Read& read);

typedef std::shared_ptr<Read> ReadPtr;

//-----------------------------------------------------------

class FastqReadArchive {
    std::vector<ReadPtr> reads_;
    std::map<std::string, ReadPtr> name_read_map_;

public:
    FastqReadArchive(std::string fastq_file_fname);

    size_t size() const;

    typedef std::vector<ReadPtr>::const_iterator read_iterator;

    read_iterator cbegin() const { return reads_.cbegin(); }

    read_iterator cend() const { return reads_.cend(); }

    ReadPtr operator[](size_t index) const;

    ReadPtr GetReadByName(std::string read_name) const;

};