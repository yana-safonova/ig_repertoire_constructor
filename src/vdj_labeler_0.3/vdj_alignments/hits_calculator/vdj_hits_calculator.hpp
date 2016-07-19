#pragma once

#include "vdj_alignments/vdj_hits_storage.hpp"
# include "read_archive.hpp"

namespace vdj_labeler {

class VDJHitsCalculator {
public:
    VDJHitsCalculator() { }

    virtual VDJHitsStoragePtr ComputeHits() const = 0;

    virtual ~VDJHitsCalculator() { }
};

} // End namespace vdj_labeler