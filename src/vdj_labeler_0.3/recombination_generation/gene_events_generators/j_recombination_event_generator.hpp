#include "shms_calculators/shm_calculator.hpp"
#include "ig_gene_recombination_event_generator.hpp"
#include "shms_calculators/shm_calculator.hpp"

namespace vdj_labeler {

class JRecombinationEventGenerator: public IgGeneRecombinationEventsGenerator {
    SHMsCalculator &shms_calculator_;
    size_t max_cleavage_;
    size_t max_palindrome_;

    recombination_utils::CleavedIgGeneAlignment GenerateCleavageEvent(
        const alignment_utils::ImmuneGeneReadAlignment &v_alignment,
        const size_t cleavage_length) const;

    void GenerateCleavageEvents(const alignment_utils::ImmuneGeneReadAlignment &v_alignment,
                                recombination_utils::IgGeneRecombinationEventStorage &v_events) const;

    recombination_utils::CleavedIgGeneAlignment GeneratePalindromicEvent(
        const alignment_utils::ImmuneGeneReadAlignment &v_alignment,
        const size_t palindrome_length) const;

    void GeneratePalindromicEvents(const alignment_utils::ImmuneGeneReadAlignment &v_alignment,
                                   recombination_utils::IgGeneRecombinationEventStorage &v_events) const;

public:
    JRecombinationEventGenerator(SHMsCalculator &shms_calculator,
                                 const size_t max_cleavage,
                                 const size_t max_palindrome) : shms_calculator_(shms_calculator),
                                                                max_cleavage_(max_cleavage),
                                                                max_palindrome_(max_palindrome) { }

    recombination_utils::IgGeneRecombinationEventStorage ComputeEvents(
        const alignment_utils::ImmuneGeneReadAlignment &gene_segment_alignment) const override;
};

} // End namespace vdj_labeler