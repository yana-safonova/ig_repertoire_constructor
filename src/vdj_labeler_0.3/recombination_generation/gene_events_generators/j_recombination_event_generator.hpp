#include <recombination_generation/gene_events_generators/shms_calculators/shm_calculator.hpp>
#include "ig_gene_recombination_event_generator.hpp"
#include "shms_calculators/shm_calculator.hpp"

class JRecombinationEventGenerator : public IgGeneRecombinationEventsGenerator {
    SHMsCalculator& shms_calculator_;
    size_t max_cleavage_;
    size_t max_palindrome_;

    CleavedIgGeneAlignment GenerateCleavageEvent(IgGeneAlignmentPtr v_alignment, size_t cleavage_length);

    void GenerateCleavageEvents(IgGeneAlignmentPtr v_alignment,
                                IgGeneRecombinationEventStoragePtr v_events);

    CleavedIgGeneAlignment GeneratePalindromicEvent(IgGeneAlignmentPtr v_alignment, size_t palindrome_length);

    void GeneratePalindromicEvents(IgGeneAlignmentPtr v_alignment,
                                   IgGeneRecombinationEventStoragePtr v_events);

public:
    JRecombinationEventGenerator(SHMsCalculator& shms_calculator,
                                 size_t max_cleavage,
                                 size_t max_palindrome) : shms_calculator_(shms_calculator),
                                                          max_cleavage_(max_cleavage),
                                                          max_palindrome_(max_palindrome) { }

    IgGeneRecombinationEventStoragePtr ComputeEvents(IgGeneAlignmentPtr gene_segment_alignment);
};