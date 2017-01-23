//
// Created by Andrew Bzikadze on 5/21/16.
//

#pragma once

#include <memory>
#include "evolutionary_edge_alignment/evolutionary_edge_alignment.hpp"

namespace shm_kmer_matrix_estimator {

class AbstractAlignmentCropper {
public:
    virtual void crop(EvolutionaryEdgeAlignment &) const = 0;
    virtual ~AbstractAlignmentCropper() { }
};
using AbstractAlignmentCropperPtr = std::unique_ptr<AbstractAlignmentCropper>;

} // End namespace shm_kmer_matrix_estimator
