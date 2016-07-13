#pragma once

#include "insertion_event_generator.hpp"

namespace vdj_labeler {

class VersatileInsertionGenerator: public InsertionEventGenerator {
public:
    recombination_utils::InsertionEventStoragePtr ComputeInsertionEvents(
        recombination_utils::CleavedIgGeneAlignment left_gene_alignment,
        recombination_utils::CleavedIgGeneAlignment right_gene_alignment);
};

} // End namespace vdj_labeler