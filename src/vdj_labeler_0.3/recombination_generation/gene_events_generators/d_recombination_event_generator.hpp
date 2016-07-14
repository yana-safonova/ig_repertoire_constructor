#pragma once

#include "ig_gene_recombination_event_generator.hpp"
#include "shms_calculators/shm_calculator.hpp"

namespace vdj_labeler {

class DRecombinationEventGenerator : public IgGeneRecombinationEventsGenerator {
    SHMsCalculator& shm_calculator_;
    size_t max_cleavage_;
    size_t max_palindrome_;

    int ComputeMinLeftBound(const alignment_utils::ImmuneGeneReadAlignmentPtr d_alignment) const;

    int ComputeMinRightBound(const alignment_utils::ImmuneGeneReadAlignmentPtr d_alignment) const;

    size_t ComputeMaxRightConsistentCleavage(const alignment_utils::ImmuneGeneReadAlignmentPtr d_alignment,
                                             const int left_event_size) const;

    void GenerateRightConsistentEvents(const alignment_utils::ImmuneGeneReadAlignmentPtr d_alignment,
                                       const int left_event_size,
                                       const recombination_utils::IgGeneRecombinationEventStoragePtr d_events) const;

public:
    DRecombinationEventGenerator(SHMsCalculator& shm_calculator,
                                 const size_t max_cleavage, const size_t max_palindrome) :
            shm_calculator_(shm_calculator),
            max_cleavage_(max_cleavage),
            max_palindrome_(max_palindrome) { }

    recombination_utils::IgGeneRecombinationEventStoragePtr ComputeEvents(
        const alignment_utils::ImmuneGeneReadAlignmentPtr gene_segment_alignment) const;
};

} // End namespace vdj_labeler
