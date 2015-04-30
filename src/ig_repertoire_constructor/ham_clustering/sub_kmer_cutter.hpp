//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "logger/logger.hpp"

#include <iostream>
#include <vector>

#include "sub_kmer_data.hpp"

namespace ham_clustering {

template <class KMerData>
class SubKMerCutterBySubPermutation {
    std::vector <unsigned> permutation_;
    unsigned from_;
    unsigned to_;
public:
    SubKMerCutterBySubPermutation(const std::vector <unsigned> & permutation, unsigned from, unsigned to)
        : permutation_(permutation), from_(from), to_(to) {
    }

    SubKMerData <KMerData> CutSubKMer(const typename KMerData::KMer & kmer, size_t kmer_index) const {
        SubKMerData <KMerData> sub_kmer;

        sub_kmer.idx = kmer_index;

        std::vector <char> chars(to_ - from_);
        for (unsigned i = from_, j = 0; i < to_; ++i, ++j) {
            chars[j] = kmer[permutation_[i]];
        }
        sub_kmer.data = KMerData::CreateSubKMer(chars);

        // Yay for NRVO!
        return sub_kmer;
    }
};

template <class KMerData>
class SubKMerCutterWithStrideFactory {
    std::vector <unsigned> permutation_;
    unsigned length_;
    unsigned tau_;
public:
    SubKMerCutterWithStrideFactory(unsigned length, unsigned stride, unsigned tau) : length_(length), tau_(tau) {
        permutation_.reserve(length);
        std::vector <bool> used(length, false);

        stride %= length;
        for (unsigned i = 0; i < length; ++i) {
            for (unsigned j = i; !used[j]; j = (j + stride) % length) {
                permutation_.push_back(j);
                used[j] = true;
            }
        }
    }

    SubKMerCutterBySubPermutation <KMerData> Create(unsigned number) const {
        unsigned from = number * length_ / (tau_ + 1);
        unsigned to = (number + 1) * length_ / (tau_ + 1);
        return SubKMerCutterBySubPermutation <KMerData>(permutation_, from, to);
    }

    unsigned GetCuttersCount() const {
        return tau_ + 1;
    }
};

template <class KMerData>
class SubKMerCutterOnlyWithMismatchesFactory {
    std::vector <unsigned> permutation_;
    unsigned length_;
    unsigned tau_;
public:
    SubKMerCutterOnlyWithMismatchesFactory(const KMerData & data, const std::vector <size_t> indices) : tau_(data.GetMaxCountMismatches()) {
        unsigned length = data.GetKmerLength();
        std::vector <int> occurred_nucls(length, 0);
        for (size_t index : indices) {
            const auto & kmer = data[index];
            for (unsigned j = 0; j < length; ++j) {
                occurred_nucls[j] |= (1 << kmer[j]);
            }
        }

        std::vector <unsigned> positions;
        positions.reserve(length);
        for (unsigned j = 0; j < length; ++j) {
            if ((occurred_nucls[j] & (occurred_nucls[j] - 1)) != 0) {
                positions.push_back(j);
            }
        }

        length_ = (unsigned)positions.size();
        std::vector <bool> used(length_, false);
        unsigned stride = tau_ + 1;

        permutation_.reserve(length_);
        for (unsigned i = 0; i < length_; ++i) {
            for (unsigned j = i; !used[j]; j = (j + stride) % length_) {
                permutation_.push_back(positions[j]);
                used[j] = true;
            }
        }
    }

    SubKMerCutterBySubPermutation <KMerData> Create(unsigned number) const {
        unsigned from = number * length_ / (tau_ + 1);
        unsigned to = (number + 1) * length_ / (tau_ + 1);
        return SubKMerCutterBySubPermutation <KMerData>(permutation_, from, to);
    }

    unsigned GetCuttersCount() const {
        return tau_ + 1;
    }

    bool HasEmptyCutters() const {
    	return length_ < GetCuttersCount();
    }
};

}
