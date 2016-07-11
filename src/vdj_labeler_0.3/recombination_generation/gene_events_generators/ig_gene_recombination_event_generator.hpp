#pragma once

#include "../../vdj_alignments/alignment_structs.hpp"
#include "../../recombination/cleaved_gene.hpp"
#include "../recombination_event_storage.hpp"

class IgGeneRecombinationEventsGenerator {
public:
    virtual IgGeneRecombinationEventStoragePtr ComputeEvents(IgGeneAlignmentPtr gene_segment_alignment) = 0;
    virtual ~IgGeneRecombinationEventsGenerator() { }
};