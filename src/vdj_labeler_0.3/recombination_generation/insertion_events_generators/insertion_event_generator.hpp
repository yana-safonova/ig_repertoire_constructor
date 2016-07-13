#pragma once

#include "recombination_utils/insertion_event_storage.hpp"
#include "recombination_utils/cleaved_gene.hpp"

namespace vdj_labeler {

class InsertionEventGenerator {
public:
    virtual recombination_utils::InsertionEventStoragePtr ComputeInsertionEvents(
        recombination_utils::CleavedIgGeneAlignment left_gene_alignment,
        recombination_utils::CleavedIgGeneAlignment right_gene_alignment) = 0;

    virtual ~InsertionEventGenerator() { }
};

} // End namespace vdj_labeler