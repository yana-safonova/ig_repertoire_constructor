#pragma once

#include <vector>
#include <algorithm>

namespace ig_repertoire_constructor {

class AlignedReadCluster {
public:
	struct ShiftedRead {
	    size_t read_number;
	    unsigned shift; // shift is non-negative. it corresponds to right-shift.
	};

	typedef std::vector <ShiftedRead>::iterator iterator;
	typedef std::vector <ShiftedRead>::const_iterator const_iterator;

	iterator begin() {
		return reads_.begin();
	}

	iterator end() {
		return reads_.end();
	}

	const_iterator begin() const {
		return reads_.begin();
	}

	const_iterator end() const {
		return reads_.end();
	}

	const ShiftedRead & operator[](size_t index) const {
	    return reads_[index];
	}

	void Add(size_t read_number, unsigned shift) {
		reads_.push_back({read_number, shift});
	}

	size_t size() const {
	    return reads_.size();
	}

private:
	vector <ShiftedRead> reads_;
};

typedef std::shared_ptr <std::vector <AlignedReadCluster> > VectorAlignedReadClusterPtr;

class AlignClusterBuilder {
public:
    AlignClusterBuilder() : minimal_shift_(0) {}

    void Add(size_t read_number, int shift) {
        minimal_shift_ = std::min(minimal_shift_, shift);
        shifted_reads_.push_back(std::make_pair(read_number, shift));
    }

    AlignedReadCluster GetCluster() const {
        AlignedReadCluster cluster;
        for (const auto & pair : shifted_reads_) {
            cluster.Add(pair.first, (unsigned)(pair.second - minimal_shift_));
        }
        return cluster;
    }

private:
    std::vector <std::pair<size_t, int> > shifted_reads_;
    int minimal_shift_;
};

}
