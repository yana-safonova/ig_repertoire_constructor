#pragma once

#include "alignment_utils/pairwise_alignment.hpp"
#include "recombination_utils/cleaved_gene.hpp"
#include "recombination_utils/recombination_event_storage.hpp"

namespace vdj_labeler {

class IgGeneRecombinationEventsGenerator {
public:
    virtual recombination_utils::IgGeneRecombinationEventStoragePtr ComputeEvents(
        alignment_utils::ImmuneGeneReadAlignmentPtr gene_segment_alignment) = 0;
    virtual ~IgGeneRecombinationEventsGenerator() { }
};

} // End namespace vdj_labeler