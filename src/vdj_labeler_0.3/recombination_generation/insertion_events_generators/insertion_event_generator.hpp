#pragma once

#include "recombination_utils/insertion_event_storage.hpp"
#include "recombination_utils/cleaved_gene.hpp"

namespace vdj_labeler {

class InsertionEventGenerator {
public:
    virtual recombination_utils::InsertionEventStoragePtr ComputeInsertionEvents(
        const recombination_utils::CleavedIgGeneAlignment &left_gene_alignment,
        const recombination_utils::CleavedIgGeneAlignment &right_gene_alignment) const = 0;

    virtual ~InsertionEventGenerator() { }
};

} // End namespace vdj_labeler