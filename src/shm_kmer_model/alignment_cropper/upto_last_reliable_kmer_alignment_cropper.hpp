//
// Created by Andrew Bzikadze on 5/21/16.
//

#pragma once

#include <string>

#include <cmath>

#include "../shm_config.hpp"
#include "abstract_alignment_cropper.hpp"

class UptoLastReliableKmerAlignmentCropper : public ns_abstract_alignment_cropper::AbstractAlignmentCropper {
private:
    const unsigned int kmer_len;
    const unsigned int hash_base;
    const unsigned int hash_max_pow;

public:
    UptoLastReliableKmerAlignmentCropper(const shm_config::alignment_cropper_params::
        upto_reliable_kmer_cropper_params& config);

    void crop(ns_gene_alignment::GermlineReadPair& alignment);
    ~UptoLastReliableKmerAlignmentCropper() { }

private:
    template <typename PairIter>
    PairIter find_correct_boarder(const PairIter&, const PairIter&);
};