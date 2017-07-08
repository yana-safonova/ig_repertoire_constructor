//
// Created by Andrew Bzikadze on 5/22/16.
//

#include <cassert>
#include "upto_last_reliable_kmer_alignment_cropper.hpp"

namespace shm_kmer_matrix_estimator {

UptoLastReliableKmerAlignmentCropper::UptoLastReliableKmerAlignmentCropper(
    const shm_kmer_matrix_estimator_config::alignment_cropper_params::upto_reliable_kmer_cropper_params &shm_config) :
    kmer_len(shm_config.kmer_len),
    hash_base(shm_config.hash_base),
    hash_max_pow(static_cast<unsigned int>(
                     std::round(
                         std::pow(hash_base, kmer_len - 1)))) {}

void UptoLastReliableKmerAlignmentCropper::crop(EvolutionaryEdgeAlignment &alignment) const {
    if (alignment.IsCropped())
        return;

    using std::pair;
    using std::make_pair;
    using std::string;

    using PairStringCIterator = pair<string::const_iterator, string::const_iterator>;
    using PairStringCRIterator = pair<string::const_reverse_iterator, string::const_reverse_iterator>;

    // Find the correct left and right boarder. After that we crop alignment up these boarders.
    auto left_boarder = find_correct_boarder<PairStringCIterator>(
        make_pair(alignment.son().cbegin(), alignment.parent().cbegin()),
        make_pair(alignment.son().cend(), alignment.parent().cend()));

    auto right_boarder = find_correct_boarder<PairStringCRIterator>(
        make_pair(alignment.son().crbegin(), alignment.parent().crbegin()),
        make_pair(alignment.son().crend(), alignment.parent().crend()));

    alignment.substract_cdr_positions(left_boarder.first - alignment.son().cbegin());

    alignment.set_son(left_boarder.first, right_boarder.first.base());
    alignment.set_parent(left_boarder.second, right_boarder.second.base());
    alignment.SetCropped();
}

template<typename PairIter>
PairIter UptoLastReliableKmerAlignmentCropper::find_correct_boarder(
    const PairIter &begin_iterators, const PairIter &end_iterators) const {
    long long hash_parent = 0, hash_son = 0, hash_mult = 1;

    auto iter_boarder(std::make_pair(begin_iterators.first + kmer_len,
                                     begin_iterators.second + kmer_len));
    // Calculate hashes of first kmer in both parent gene and son.
    while (iter_boarder.first != begin_iterators.first) {
        --iter_boarder.first;
        --iter_boarder.second;

        hash_parent += (*iter_boarder.first) * hash_mult;
        hash_son += (*iter_boarder.second) * hash_mult;

        hash_mult *= hash_base;
    }
    // If the hashes are equal — no need to crop anything.
    if (hash_parent == hash_son)
        return begin_iterators;

    assert(iter_boarder.first == begin_iterators.first);
    assert(iter_boarder.second == begin_iterators.second);

    // Otherwise, find a kmer with equal hashes in the parent gene sequence and the son.
    while (iter_boarder.first + kmer_len != end_iterators.first &&
        hash_parent != hash_son) {
        hash_parent -= (*iter_boarder.first) * hash_max_pow;
        hash_son -= (*iter_boarder.second) * hash_max_pow;

        hash_parent *= hash_base;
        hash_son *= hash_base;

        hash_parent += *(iter_boarder.first + kmer_len);
        hash_son += *(iter_boarder.second + kmer_len);

        ++iter_boarder.first;
        ++iter_boarder.second;
    }
    return iter_boarder;
}

} // End namespace shm_kmer_matrix_estimator
