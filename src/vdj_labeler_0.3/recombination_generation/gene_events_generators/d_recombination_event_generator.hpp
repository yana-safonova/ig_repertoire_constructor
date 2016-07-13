#pragma once

#include "ig_gene_recombination_event_generator.hpp"
#include "shms_calculators/shm_calculator.hpp"

namespace vdj_labeler {

class DRecombinationEventGenerator : public IgGeneRecombinationEventsGenerator {
    SHMsCalculator& shm_calculator_;
    size_t max_cleavage_;
    size_t max_palindrome_;

    int ComputeMinLeftBound(alignment_utils::ImmuneGeneReadAlignmentPtr d_alignment);

    int ComputeMinRightBound(alignment_utils::ImmuneGeneReadAlignmentPtr d_alignment);

    size_t ComputeMaxRightConsistentCleavage(alignment_utils::ImmuneGeneReadAlignmentPtr d_alignment,
                                          int left_event_size);

    void GenerateRightConsistentEvents(alignment_utils::ImmuneGeneReadAlignmentPtr d_alignment,
                                       int left_event_size,
                                       recombination_utils::IgGeneRecombinationEventStoragePtr d_events);

public:
    DRecombinationEventGenerator(SHMsCalculator& shm_calculator,
                                 size_t max_cleavage, size_t max_palindrome) :
            shm_calculator_(shm_calculator),
            max_cleavage_(max_cleavage),
            max_palindrome_(max_palindrome) { }

    recombination_utils::IgGeneRecombinationEventStoragePtr ComputeEvents(
        alignment_utils::ImmuneGeneReadAlignmentPtr gene_segment_alignment);
};

} // End namespace vdj_labeler
