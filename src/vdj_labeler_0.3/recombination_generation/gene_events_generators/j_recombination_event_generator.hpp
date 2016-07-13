#include "shms_calculators/shm_calculator.hpp"
#include "ig_gene_recombination_event_generator.hpp"
#include "shms_calculators/shm_calculator.hpp"

namespace vdj_labeler {

class JRecombinationEventGenerator: public IgGeneRecombinationEventsGenerator {
    SHMsCalculator &shms_calculator_;
    size_t max_cleavage_;
    size_t max_palindrome_;

    recombination_utils::CleavedIgGeneAlignment GenerateCleavageEvent(
        alignment_utils::ImmuneGeneReadAlignmentPtr v_alignment,
        size_t cleavage_length);

    void GenerateCleavageEvents(alignment_utils::ImmuneGeneReadAlignmentPtr v_alignment,
                                recombination_utils::IgGeneRecombinationEventStoragePtr v_events);

    recombination_utils::CleavedIgGeneAlignment GeneratePalindromicEvent(
        alignment_utils::ImmuneGeneReadAlignmentPtr v_alignment,
        size_t palindrome_length);

    void GeneratePalindromicEvents(alignment_utils::ImmuneGeneReadAlignmentPtr v_alignment,
                                   recombination_utils::IgGeneRecombinationEventStoragePtr v_events);

public:
    JRecombinationEventGenerator(SHMsCalculator &shms_calculator,
                                 size_t max_cleavage,
                                 size_t max_palindrome) : shms_calculator_(shms_calculator),
                                                          max_cleavage_(max_cleavage),
                                                          max_palindrome_(max_palindrome) { }

    recombination_utils::IgGeneRecombinationEventStoragePtr ComputeEvents(
        alignment_utils::ImmuneGeneReadAlignmentPtr gene_segment_alignment);
};

} // End namespace vdj_labeler