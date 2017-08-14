//
// Created by Andrew Bzikadze on 5/21/16.
//

#pragma once

#include <string>

#include <cmath>

#include "shm_kmer_matrix_estimator_config.hpp"
#include "abstract_alignment_cropper.hpp"

namespace shm_kmer_matrix_estimator {

class UptoLastReliableKmerAlignmentCropper: public AbstractAlignmentCropper {
private:
    const unsigned int kmer_len;
    const unsigned int hash_base;
    const unsigned int hash_max_pow;

public:
    explicit UptoLastReliableKmerAlignmentCropper(const shm_kmer_matrix_estimator_config::alignment_cropper_params::
                                                  upto_reliable_kmer_cropper_params &config);

    void crop(EvolutionaryEdgeAlignment &alignment) const override;
    virtual ~UptoLastReliableKmerAlignmentCropper() { }

private:
    template<typename PairIter>
    PairIter find_correct_boarder(const PairIter &, const PairIter &) const;
};

} // End namespace shm_kmer_matrix_estimator