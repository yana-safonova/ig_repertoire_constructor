#pragma once

#include "vdj_alignments/vdj_hits.hpp"

namespace vdj_labeler {

template<class Recombination, class RecombinationStoragePtr>
class AbstractRecombinationGenerator {
public:
    virtual RecombinationStoragePtr ComputeRecombinations(VDJHitsPtr vdj_hits) = 0;

    virtual ~AbstractRecombinationGenerator() { }
};

} // End namespace vdj_labeler