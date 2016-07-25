#pragma once

#include "recombination_utils/insertion_event_storage.hpp"
#include "recombination_utils/cleaved_gene.hpp"

namespace vdj_labeler {

class InsertionEventGenerator {
public:
    InsertionEventGenerator() = default;

    InsertionEventGenerator(const InsertionEventGenerator &)           = delete;
    InsertionEventGenerator& operator=(const InsertionEventGenerator&) = delete;
    InsertionEventGenerator(InsertionEventGenerator &&)                = delete;
    InsertionEventGenerator& operator=(InsertionEventGenerator&&)      = delete;

    virtual recombination_utils::InsertionEventStorage ComputeInsertionEvents(
        const recombination_utils::CleavedIgGeneAlignment &left_gene_alignment,
        const recombination_utils::CleavedIgGeneAlignment &right_gene_alignment) const = 0;

    virtual ~InsertionEventGenerator() { }
};

} // End namespace vdj_labeler