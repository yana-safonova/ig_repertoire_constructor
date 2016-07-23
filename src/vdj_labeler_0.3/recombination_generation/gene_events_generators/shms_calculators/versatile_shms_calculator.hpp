#include "shm_calculator.hpp"

namespace vdj_labeler {

class VersatileGeneSHMsCalculator : public SHMsCalculator {
    SHMsCalculator& left_shms_calculator_;
    SHMsCalculator& right_shms_calculator_;

public:
    VersatileGeneSHMsCalculator(SHMsCalculator &left_shms_calculator,
                                SHMsCalculator &right_shms_calculator) :
            left_shms_calculator_(left_shms_calculator),
            right_shms_calculator_(right_shms_calculator)
    { }

    int ComputeNumberSHMs(const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                          const int left_cleavage_length,
                          const int right_cleavage_length) const;

    int ComputeNumberSHMsForLeftEvent(const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                      const int left_cleavage_length) const;

    int ComputeNumberSHMsForRightEvent(const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                       const int right_cleavage_length) const;
};

} // End namespace vdj_labeler;