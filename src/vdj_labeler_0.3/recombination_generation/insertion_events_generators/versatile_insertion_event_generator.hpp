#pragma once

#include "insertion_event_generator.hpp"

namespace vdj_labeler {

class VersatileInsertionGenerator: public InsertionEventGenerator {
public:
    recombination_utils::InsertionEventStorage ComputeInsertionEvents(
        const recombination_utils::CleavedIgGeneAlignment &left_gene_alignment,
        const recombination_utils::CleavedIgGeneAlignment &right_gene_alignment) const override;
};

} // End namespace vdj_labeler