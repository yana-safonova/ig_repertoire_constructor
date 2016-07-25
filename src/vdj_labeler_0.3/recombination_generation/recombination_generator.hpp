#pragma once

#include "vdj_alignments/vdj_hits.hpp"

namespace vdj_labeler {

template<class Recombination, class RecombinationStorage>
class AbstractRecombinationGenerator {
public:
    AbstractRecombinationGenerator() = default;

    AbstractRecombinationGenerator(const AbstractRecombinationGenerator &)           = delete;
    AbstractRecombinationGenerator& operator=(const AbstractRecombinationGenerator&) = delete;
    AbstractRecombinationGenerator(AbstractRecombinationGenerator &&)                = delete;
    AbstractRecombinationGenerator& operator=(AbstractRecombinationGenerator&&)      = delete;

    virtual RecombinationStorage ComputeRecombinations(const VDJHits &vdj_hits) = 0;

    virtual ~AbstractRecombinationGenerator() { }
};

} // End namespace vdj_labeler