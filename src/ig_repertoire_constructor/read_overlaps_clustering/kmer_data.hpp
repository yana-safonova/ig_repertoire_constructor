#pragma once

#include "spliced_read.hpp"
#include "ham_clustering/ham_clusterer.hpp"

namespace ig_repertoire_constructor {

#define MAX_KMER_LENGTH 500

class SplicedKMerData: public ham_clustering::RtSeqKMerData<MAX_KMER_LENGTH, 0, SplicedRead> {
public:
    explicit SplicedKMerData(const std::vector<SplicedRead> & splicedReads, unsigned max_count_mismatches)
            : spliced_reads_(splicedReads), max_count_mismatches_(max_count_mismatches) {
    }

    virtual KMer operator[](size_t index) const {
        return spliced_reads_[index];
    }

    virtual size_t size() const {
        return spliced_reads_.size();
    }

    virtual unsigned GetKmerLength() const {
        if (spliced_reads_.empty()) {
            return 0;
	}
        return (unsigned)spliced_reads_[0].size();
    }

    virtual unsigned GetMaxCountMismatches() const {
	return max_count_mismatches_;
    }

    virtual ~SplicedKMerData() {
    }

private:
    const std::vector<SplicedRead> & spliced_reads_;
    const unsigned max_count_mismatches_;
};

}
