#pragma once

#include "alignment_utils/pairwise_alignment.hpp"
#include "recombination_utils/cleaved_gene.hpp"
#include "recombination_utils/recombination_event_storage.hpp"

namespace vdj_labeler {

class IgGeneRecombinationEventsGenerator {
public:
    IgGeneRecombinationEventsGenerator() = default;

    IgGeneRecombinationEventsGenerator(const IgGeneRecombinationEventsGenerator &)           = delete;
    IgGeneRecombinationEventsGenerator& operator=(const IgGeneRecombinationEventsGenerator&) = delete;
    IgGeneRecombinationEventsGenerator(IgGeneRecombinationEventsGenerator &&)                = delete;
    IgGeneRecombinationEventsGenerator& operator=(IgGeneRecombinationEventsGenerator&&)      = delete;

    virtual recombination_utils::IgGeneRecombinationEventStorage ComputeEvents(
        const alignment_utils::ImmuneGeneReadAlignment &gene_segment_alignment) const = 0;
    virtual ~IgGeneRecombinationEventsGenerator() { }
};

} // End namespace vdj_labeler