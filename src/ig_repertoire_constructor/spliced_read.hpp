#pragma once

#include "io/single_read.hpp"
#include <fstream>

namespace ig_repertoire_constructor {

class SplicedRead {
public:
    SplicedRead(std::shared_ptr <ReadArchive> read_archive_ptr, size_t read_number, unsigned from, unsigned to)
        : read_archive_ptr_(read_archive_ptr), read_number_(read_number), from_(from), to_(to) {
    }

    unsigned size() const {
        return to_ - from_;
    }

    bool HasCharWithIndex(int index) const {
        return (index + (int)from_ >= 0) && (unsigned(index + (int)from_) < GetRead().size());
    }

    char operator[](long int index) const {
        return GetRead()[index + from_];
    }

    char GetQuality(long int index) const {
        return GetRead().GetQualityString()[index + from_];
    }

    std::string GetReadName() const {
        return GetRead().name();
    }

    size_t GetReadNumber() const {
        return read_number_;
    }

    unsigned GetFrom() const {
        return from_;
    }

    unsigned GetTo() const {
        return to_;
    }

    size_t AbsoletePosition(size_t index) const {
        return GetFrom() + index;
    }

    size_t ReadLength() const {
        return GetRead().size();
    }

    string ReadName() const {
        return GetRead().name();
    }

    string FullSequence() const {
        return GetRead().sequence().str();
    }

private:
    const io::SingleRead & GetRead() const {
        return read_archive_ptr_->operator[](read_number_);
    }

    std::shared_ptr <ReadArchive> read_archive_ptr_;
    size_t read_number_;
    unsigned from_;
    unsigned to_;
};

typedef std::shared_ptr <std::vector <std::vector<SplicedRead> > > VectorSplicedReadClusterPtr;

inline VectorSplicedReadClusterPtr Deserialize(const std::string & fname, std::shared_ptr <ReadArchive> read_archive_ptr) {
    MMappedReader reader(fname);
    VERIFY(reader.good());

    VectorSplicedReadClusterPtr clusters_ptr(new std::vector <std::vector <SplicedRead> >());
    while(reader.good()) {
        size_t sz;
        reader.read((char*)&sz, sizeof(sz));
        clusters_ptr->push_back(std::vector<SplicedRead>());

        std::vector<SplicedRead> &current_cluster = clusters_ptr->back();
        current_cluster.reserve(sz);
        for (size_t i = 0; i < sz; ++i) {
            size_t read_number = 0;
            unsigned from = 0, to = 0;
            reader.read((char*)&read_number, sizeof(read_number));
            reader.read((char*)&from, sizeof(from));
            reader.read((char*)&to, sizeof(to));
            current_cluster.emplace_back(read_archive_ptr, read_number, from, to);
        }
    }

    return clusters_ptr;
}

inline void Serialize(const std::string & fname, VectorSplicedReadClusterPtr clusters_ptr) {
    std::ofstream writer(fname, std::ios::out | std::ios::binary);
    VERIFY(writer.good());
    for (const auto & cluster : *clusters_ptr) {
        size_t sz = cluster.size();
        writer.write((char*)&sz, sizeof(sz));

        for (size_t i = 0; i < sz; ++i) {
            size_t read_number = cluster[i].GetReadNumber();
            unsigned from = cluster[i].GetFrom();
            unsigned to = cluster[i].GetTo();
            writer.write((char*)&read_number, sizeof(read_number));
            writer.write((char*)&from, sizeof(from));
            writer.write((char*)&to, sizeof(to));
        }
    }
    VERIFY(!writer.fail());
    writer.close();
}

inline void SerializeHumanReadable(const std::string & fname, VectorSplicedReadClusterPtr clusters_ptr) {
    std::ofstream writer(fname, std::ios::out | std::ios::binary);
    VERIFY(writer.good());
    for (const auto & cluster : *clusters_ptr) {
        for (size_t i = 0; i < cluster.size(); ++i) {
            unsigned from = cluster[i].GetFrom();
            unsigned to = cluster[i].GetTo();
            writer << cluster[i].GetReadName() << " from=" << from << " to=" << to << std::endl;
        }
        writer << "__END__OF__CLUSTER__" << std::endl;
    }
    VERIFY(!writer.fail());
    writer.close();
}

int GetLeftmostPosition(std::vector<SplicedRead> const &spliced_reads,
                        std::vector <size_t> const &subcluster) {
    for (int pos = -1; ; --pos) {
        bool has_reads = false;
        size_t max_i = (subcluster.empty()) ? spliced_reads.size() : subcluster.size();
        for (size_t i = 0; i != max_i; ++i) {
            size_t spliced_read_ind = (subcluster.empty()) ? i : subcluster[i];
            SplicedRead const &curr_read = spliced_reads[spliced_read_ind];
            if (curr_read.HasCharWithIndex(pos)) {
                has_reads = true;
                break;
            }
        }
        if (!has_reads) {
            return pos + 1;
        }
    }
}

int GetRightmostPosition(std::vector<SplicedRead> const &spliced_reads,
                        std::vector <size_t> const &subcluster) {
    for (size_t pos = spliced_reads.size(); ; ++pos) {
        bool has_reads = false;
        size_t max_i = (subcluster.empty()) ? spliced_reads.size() : subcluster.size();
        for (size_t i = 0; i != max_i; ++i) {
            size_t spliced_read_ind = (subcluster.empty()) ? i : subcluster[i];
            SplicedRead const &curr_read = spliced_reads[spliced_read_ind];
            if (curr_read.HasCharWithIndex((int)pos)) {
                has_reads = true;
                break;
            }
        }
        if (!has_reads) {
            return (int)(pos - 1);
        }
    }
}

}
