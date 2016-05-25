//
// Created by Andrew Bzikadze on 5/22/16.
//

#include "upto_last_reliable_kmer_alignment_cropper.hpp"

UptoLastReliableKmerAlignmentCropper::UptoLastReliableKmerAlignmentCropper(
    const shm_config::alignment_cropper_params::upto_reliable_kmer_cropper_params &config) :
    kmer_len(config.kmer_len),
    hash_base(config.hash_base),
    hash_max_pow(static_cast<unsigned int>(
                     std::round(
                         std::pow(hash_base, kmer_len - 1))))
    { }

void UptoLastReliableKmerAlignmentCropper::crop(ns_gene_alignment::ReadGermlineAlignment &alignment) const {
    using std::pair;
    using std::make_pair;
    using std::string;

    using PairStringCIterator = pair<string::const_iterator, string::const_iterator>;
    using PairStringCRIterator = pair<string::const_reverse_iterator, string::const_reverse_iterator>;

    auto left_boarder = find_correct_boarder<PairStringCIterator>(
        make_pair(alignment.read().cbegin(), alignment.germline().cbegin()),
        make_pair(alignment.read().cend(),   alignment.germline().cend()));

    auto right_boarder = find_correct_boarder<PairStringCRIterator>(
        make_pair(alignment.read().crbegin(), alignment.germline().crbegin()),
        make_pair(alignment.read().crend(),   alignment.germline().crend()));

    alignment.set_read(left_boarder.first, right_boarder.first.base());
    alignment.set_germline(left_boarder.second, right_boarder.second.base());
}

template<typename PairIter>
PairIter UptoLastReliableKmerAlignmentCropper::find_correct_boarder(
    const PairIter &begin_iterators, const PairIter &end_iterators) const {
    long long hash_germline = 0, hash_read = 0, hash_mult = 1;

    auto iter_boarder(std::make_pair(begin_iterators.first + kmer_len,
                                     begin_iterators.second + kmer_len));
    while (iter_boarder.first != begin_iterators.first) {
        --iter_boarder.first;
        --iter_boarder.second;

        hash_germline += (*iter_boarder.first) * hash_mult;
        hash_read += (*iter_boarder.second) * hash_mult;

        hash_mult *= hash_base;
    }
    if (hash_germline == hash_read)
        return begin_iterators;

    assert(iter_boarder.first == begin_iterators.first);
    assert(iter_boarder.second == begin_iterators.second);

    while (iter_boarder.first + kmer_len != end_iterators.first &&
        hash_germline != hash_read) {
        hash_germline -= (*iter_boarder.first) * hash_max_pow;
        hash_read -= (*iter_boarder.second) * hash_max_pow;

        hash_germline *= hash_base;
        hash_read *= hash_base;

        hash_germline += *(iter_boarder.first + kmer_len);
        hash_read += *(iter_boarder.second + kmer_len);

        ++iter_boarder.first;
        ++iter_boarder.second;
    }
    return iter_boarder;
}