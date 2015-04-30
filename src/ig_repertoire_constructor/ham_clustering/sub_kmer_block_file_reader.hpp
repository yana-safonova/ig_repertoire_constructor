#pragma once

#include "io/mmapped_reader.hpp"
#include <string>
#include <vector>

#include "sub_kmer_data.hpp"

namespace ham_clustering {

template <class KMerData>
class SubKMerBlockFileReader {
    MMappedReader ifs_;

 public:
    SubKMerBlockFileReader(const std::string &fname, bool unlink = false)
            : ifs_(fname, unlink) { }

    bool GetBlock(std::vector<size_t> &block) {
        block.clear();
#if 0
        block.shrink_to_fit();
#else
        std::vector<size_t>().swap(block);
#endif

        if (!ifs_.good())
            return false;

        size_t sz;
        ifs_.read((char*)&sz, sizeof(sz));
        block.resize(sz);
        for (size_t i = 0; i < sz; ++i) {
            SubKMerData <KMerData> s;
            binary_read(ifs_, s);
            block[i] = s.idx;
        }

        return true;
    }
};

}
