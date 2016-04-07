#pragma once

#include <unordered_map>
#include <seqan/sequence.h>

namespace core {

    struct Read {
        std::string name;
        seqan::Dna5String seq;
        size_t id;

        Read(std::string new_name,
             seqan::Dna5String new_seq,
             size_t new_id) :
                name(new_name),
                seq(new_seq),
                id(new_id) { }

        Read() : name(), seq() { }

        size_t length() const { return seqan::length(seq); }
    };

    std::ostream &operator<<(std::ostream &out, const Read &read);

    typedef std::shared_ptr <Read> ReadPtr;

//-----------------------------------------------------------

    class ReadArchive {
        std::vector <Read> reads_;
        std::unordered_map <std::string, size_t> name_index_map_;

    public:
        ReadArchive(std::string fastq_file_fname);

        size_t size() const;

        typedef std::vector<Read>::const_iterator read_iterator;

        read_iterator cbegin() const { return reads_.cbegin(); }

        read_iterator cend() const { return reads_.cend(); }

        const Read& operator[](size_t index) const;

        const Read& GetReadByName(std::string read_name) const;

        size_t GetIndexByReadName(std::string read_name) const;

        void FixSpacesInHeaders();
    };
}