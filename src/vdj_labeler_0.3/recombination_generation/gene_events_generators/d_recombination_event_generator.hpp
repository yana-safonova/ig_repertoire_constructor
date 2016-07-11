#pragma once

#include "ig_gene_recombination_event_generator.hpp"
#include "shms_calculators/shm_calculator.hpp"

class DRecombinationEventGenerator : public IgGeneRecombinationEventsGenerator {
    SHMsCalculator& shm_calculator_;
    size_t max_cleavage_;
    size_t max_palindrome_;

    int ComputeMinLeftBound(IgGeneAlignmentPtr d_alignment);

    int ComputeMinRightBound(IgGeneAlignmentPtr d_alignment);

    size_t ComputeMaxRightConsistentCleavage(IgGeneAlignmentPtr d_alignment,
                                          int left_event_size);

    void GenerateRightConsistentEvents(IgGeneAlignmentPtr d_alignment,
                                       int left_event_size,
                                       IgGeneRecombinationEventStoragePtr d_events);

public:
    DRecombinationEventGenerator(SHMsCalculator& shm_calculator,
                                 size_t max_cleavage, size_t max_palindrome) :
            shm_calculator_(shm_calculator),
            max_cleavage_(max_cleavage),
            max_palindrome_(max_palindrome) { }

    IgGeneRecombinationEventStoragePtr ComputeEvents(IgGeneAlignmentPtr gene_segment_alignment);
};