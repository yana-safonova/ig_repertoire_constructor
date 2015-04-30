#pragma once

#include "logger/logger.hpp"
#include "io/single_read.hpp"

// TODO all these includes for one method: LoadReadArchive()
#include "standard.hpp"
#include "io/library.hpp"
#include "io/io_helper.hpp"
#include "ig_config.hpp"

namespace ig_repertoire_constructor {

class ReadArchive {
public:
    typedef std::vector <io::SingleRead>::iterator iterator;
    typedef std::vector <io::SingleRead>::const_iterator const_iterator;

    io::SingleRead& operator[](size_t number);
    const io::SingleRead& operator[](size_t number) const;
    size_t GetReadNumberByReadName(const std::string & read_name) const;

    size_t size() const;

    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;

    template <typename Stream>
    void LoadFromStream(Stream & stream);

private:
    void AddNewRead(const io::SingleRead & read);

    std::map <std::string, size_t> read_number_by_name_;
    std::vector <io::SingleRead> reads_;
};

typedef std::shared_ptr <ReadArchive> ReadArchivePtr;

template <typename Stream>
inline void ReadArchive::LoadFromStream(Stream & stream) {
    reads_.resize(0);
    io::SingleRead read;
    int skipped_reads = 0;
    while (!stream.eof()){
        stream >> read;
        if (read.IsValid()) {
            AddNewRead(read);
        } else {
            skipped_reads++;
        }
    }
    if (skipped_reads != 0) {
        INFO(skipped_reads << " reads are invalid. They were been skipped.");
    }
    INFO("Total count of reads: " << reads_.size());
}

inline std::shared_ptr <ReadArchive> LoadReadArchive() {
    const io::DataSet <> &dataset = ig_cfg::get().io.dataset;
    io::ReadStreamList <io::SingleRead> readers;
    for (auto I = dataset.reads_begin(), E = dataset.reads_end(); I != E; ++I) {
        readers.push_back(io::EasyStream(*I, false));
    }
    auto stream_ptr = MultifileWrap(readers);
    auto read_archive = std::make_shared <ig_repertoire_constructor::ReadArchive>();
    read_archive->LoadFromStream(*stream_ptr);
    return read_archive;
}

}
