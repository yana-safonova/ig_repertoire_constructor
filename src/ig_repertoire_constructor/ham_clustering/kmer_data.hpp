#pragma once

#include "sequence/rtseq.hpp"

namespace ham_clustering {

template <size_t MAX_KMER_LENGTH, size_t MIN_TAU, typename KMER_TYPE_NAME = RuntimeSeq<MAX_KMER_LENGTH> >
class RtSeqKMerData {
public:
    typedef KMER_TYPE_NAME KMer;
    typedef RuntimeSeq<(MAX_KMER_LENGTH + MIN_TAU) / (MIN_TAU + 1)> SubKMer;

    template <class Reader>
    static void binary_read_subkmer(Reader &is, SubKMer &sub_kmer) {
        typename SubKMer::DataType seq_data[SubKMer::DataSize];

        short seq_size;
        is.read((char*)&seq_size, sizeof(seq_size));
        is.read((char*)seq_data, sizeof(seq_data));

        sub_kmer = SubKMer((size_t)seq_size, seq_data);
    }

    template <class Writer>
    static void binary_write_subkmer(Writer &os, const SubKMer &sub_kmer) {
        short seq_size = (short)sub_kmer.size();
        os.write((char*)&seq_size, sizeof(seq_size));
        os.write((char*)sub_kmer.data(), SubKMer::TotalBytes);
    }

    static SubKMer CreateSubKMer(const std::vector<char> & chars) {
        return SubKMer(chars.size(), chars);
    }

    virtual KMer operator[](size_t index) const = 0;
    virtual size_t size() const = 0;

    virtual unsigned GetKmerLength() const = 0;
    virtual unsigned GetMaxCountMismatches() const = 0;

    virtual ~RtSeqKMerData() {
    }
};

}
