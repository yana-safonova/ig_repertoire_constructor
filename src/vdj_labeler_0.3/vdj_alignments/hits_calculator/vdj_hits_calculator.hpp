#pragma once

#include "vdj_alignemnts/vdj_hits_storage.hpp"
# include "read_archive.hpp"

namespace vdj_labeler {

class VDJHitsCalculator {
protected:
    const FastqReadArchive &read_archive_;

public:
    VDJHitsCalculator(const FastqReadArchive &read_archive) :
        read_archive_(read_archive) { }

    virtual VDJHitsStoragePtr ComputeHits() = 0;

    virtual ~VDJHitsCalculator() { }
};

} // End namespace vdj_labeler