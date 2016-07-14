#include "versatile_shms_calculator.hpp"

using namespace vdj_labeler;

int VersatileGeneSHMsCalculator::ComputeNumberSHMs(const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
                                                   const int left_cleavage_length,
                                                   const int right_cleavage_length) const {
    return int(gene_alignment->NumberSHMs()) +
            left_shms_calculator_.ComputeNumberSHMs(gene_alignment, left_cleavage_length, right_cleavage_length) +
            right_shms_calculator_.ComputeNumberSHMs(gene_alignment, left_cleavage_length, right_cleavage_length);
}

int VersatileGeneSHMsCalculator::ComputeNumberSHMsForLeftEvent(
        const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
        const int left_cleavage_length) const
{
    return left_shms_calculator_.ComputeNumberSHMsForLeftEvent(gene_alignment, left_cleavage_length);
}

int VersatileGeneSHMsCalculator::ComputeNumberSHMsForRightEvent(
        const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
        const int right_cleavage_length) const
{
    return left_shms_calculator_.ComputeNumberSHMsForRightEvent(gene_alignment, right_cleavage_length);
}