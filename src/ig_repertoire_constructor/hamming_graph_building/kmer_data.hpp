#pragma once

#include "ham_clustering/kmer_data.hpp"
#include <algorithm>
#include <vector>

namespace ig_repertoire_constructor {

#define MAX_KMER_LENGTH 500
#define MAX_MISMATCHES_ON_OVERLAP 3

class KMerInfo;

class KMerData: public ham_clustering::RtSeqKMerData<MAX_KMER_LENGTH, MAX_MISMATCHES_ON_OVERLAP> {
public:
    KMerData(const ReadArchive & read_archive, unsigned kmer_length);

    const KMerInfo & GetInfo(size_t index) const;
    std::pair<size_t, size_t> GetRange(size_t read_index) const;

    virtual KMer operator[](size_t index) const;
    virtual size_t size() const;

    virtual unsigned GetKmerLength() const {
        return kmer_length_;
    }

    virtual unsigned GetMaxCountMismatches() const {
        return MAX_MISMATCHES_ON_OVERLAP;
    }

    virtual ~KMerData() {
    }
private:
    unsigned kmer_length_;
    std::vector<KMerInfo> kmer_infos_;
    std::vector<std::pair<size_t, size_t> > ranges_;
};

class KMerInfo {
public:
    KMerInfo(size_t read_number, unsigned offset, const KMerData::KMer & kmer)
            : read_number(read_number), offset(offset), kmer(kmer) {
    }

    size_t read_number;
    unsigned offset;
    KMerData::KMer kmer;
};

KMerData::KMerData(const ReadArchive & read_archive, unsigned kmer_length) : kmer_length_(kmer_length) {
    ranges_.resize(read_archive.size());
    for (size_t index = 0; index < read_archive.size(); ++index) {
        ranges_[index].first = ranges_[index].second = kmer_infos_.size();

        Sequence sequence = read_archive[index].sequence();
        if (sequence.size() < kmer_length_) {
            continue;
        }
        KMerData::KMer kmer = sequence.start<KMerData::KMer>(kmer_length_);
        kmer_infos_.emplace_back(index, 0, kmer);
        for (unsigned offset = 1; offset + kmer_length_ <= sequence.size(); ++offset) {
            kmer <<= sequence[offset + kmer_length_ - 1];
            kmer_infos_.emplace_back(index, offset, kmer);
        }

        ranges_[index].second = kmer_infos_.size();
    }
}

const KMerInfo & KMerData::GetInfo(size_t index) const {
    return kmer_infos_[index];
}

std::pair<size_t, size_t> KMerData::GetRange(size_t read_index) const {
    return ranges_[read_index];
}

KMerData::KMer KMerData::operator[](size_t index) const {
    return kmer_infos_[index].kmer;
}

size_t KMerData::size() const {
    return kmer_infos_.size();
}

}
