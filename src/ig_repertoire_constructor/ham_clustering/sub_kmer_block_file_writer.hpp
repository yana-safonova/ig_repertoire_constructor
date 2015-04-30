#pragma once

#include <fstream>
#include <string>
#include <vector>
#include "verify.hpp"

#include "sub_kmer_data.hpp"
#include "sub_kmer_cutter.hpp"

namespace ham_clustering {

template <class KMerData>
class SubKMerBlockFileWriter {
public:
    SubKMerBlockFileWriter(const std::string & fname, const KMerData & data)
        : fname_(fname), data_(data), is_opened_(false), is_closed_(false), count_written_blocks_(0) {
    }

    size_t GetCountWrittenBlocks() const {
        return count_written_blocks_;
    }

    template <class SubKMerCutter>
    void Write(const SubKMerCutter &cutter, const std::vector<size_t> *block) {
        EnsureOpenStream();

        size_t sz = (block == NULL ? data_.size() : block->size());
        ofs_.write((char*)&sz, sizeof(sz));
        for (size_t i = 0, e = sz; i != e; ++i) {
            size_t idx = (block == NULL ? i : (*block)[i]);
            SubKMerData <KMerData> s = cutter.CutSubKMer(data_[idx], idx);
            binary_write(ofs_, s);
        }
        count_written_blocks_++;
    }

    void Close() {
    	if (is_closed_) {
    		return;
    	}
        if (is_opened_) {
            VERIFY(!ofs_.fail());
            ofs_.flush();
            ofs_.close();
            is_closed_ = true;
        }
    }

    ~SubKMerBlockFileWriter() {
        Close();
    }

private:

    void EnsureOpenStream() {
    	VERIFY_MSG(!is_closed_, "File '" << fname_ << "' is already closed!");

        if (!is_opened_) {
            is_opened_ = true;
            ofs_.open(fname_, std::ios::out | std::ios::binary);
            VERIFY(ofs_.good());
        }
    }

    std::string fname_;
    const KMerData & data_;
    bool is_opened_;
    bool is_closed_;
    size_t count_written_blocks_;
    std::ofstream ofs_;
};

}
