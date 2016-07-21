#include <cmath>

#include "hc_model_based_recombination_calculator.hpp"

namespace vdj_labeler {

double HCModelBasedRecombinationCalculator::ComputeAssemblyProbability(
        const recombination_utils::HCRecombination &recombination) const
{
    /** Probabilities of Gene Ids are computed */
    double V_gene_prob_id = model_.GetProbabilityByGenId(germline_utils::SegmentType::VariableSegment,
                                                         recombination.V().GeneId());
    V_gene_prob_id = log(V_gene_prob_id);

    double D_gene_prob_id = model_.GetProbabilityByGenId(germline_utils::SegmentType::DiversitySegment,
                                                         recombination.D().GeneId());
    D_gene_prob_id = log(D_gene_prob_id);

    double J_gene_prob_id = model_.GetProbabilityByGenId(germline_utils::SegmentType::JoinSegment,
                                                         recombination.J().GeneId());
    J_gene_prob_id = log(J_gene_prob_id);

    double id_probs = V_gene_prob_id + D_gene_prob_id + J_gene_prob_id;

    /** Deletion probabilities */
    double V_gene_prob_del = model_.GetVPalindromeDeletionModel().GetDeletionProbability(
        recombination.V().GeneId(),
        recombination.V().RightCleavageLength());
    V_gene_prob_del = log(V_gene_prob_del);

    double DLeft_gene_prob_del = model_.GetDLeftPalindromeDeletionModel().GetDeletionProbability(
        recombination.D().GeneId(),
        recombination.D().LeftCleavageLength());
    double DRight_gene_prob_del = model_.GetDRightPalindromeDeletionModel().GetDeletionProbability(
        recombination.D().GeneId(),
        recombination.D().RightCleavageLength());
    DLeft_gene_prob_del = log(DLeft_gene_prob_del);
    DRight_gene_prob_del = log(DRight_gene_prob_del);

    double J_gene_prob_del = model_.GetJPalindromeDeletionModel().GetDeletionProbability(
        recombination.J().GeneId(),
        recombination.J().LeftCleavageLength());
    J_gene_prob_del = log(J_gene_prob_del);

    double del_probs = V_gene_prob_del + DLeft_gene_prob_del +
        DRight_gene_prob_del + J_gene_prob_del;

    /** Insertion probabilities */
    double VD_non_genomic_ins_prob = model_.GetVDNongenomicInsertionModel().GetInsertionProbabilityByLength(
        recombination.VDInsertion().length());
    VD_non_genomic_ins_prob = log(VD_non_genomic_ins_prob);
    double DJ_non_genomic_ins_prob = model_.GetDJNongenomicInsertionModel().GetInsertionProbabilityByLength(
        recombination.DJInsertion().length());
    DJ_non_genomic_ins_prob = log(DJ_non_genomic_ins_prob);

    for (size_t i = 0; i < recombination.VDInsertion().length(); ++i)
        VD_non_genomic_ins_prob += log(model_.GetVDNongenomicInsertionModel().GetTransitionProbability(
            {recombination.Read().seq[i + recombination.VDInsertion().StartPosition() - 1],
             recombination.Read().seq[i + recombination.VDInsertion().StartPosition()]}));

    for (size_t i = 0; i < recombination.DJInsertion().length(); ++i)
        DJ_non_genomic_ins_prob += log(model_.GetDJNongenomicInsertionModel().GetTransitionProbability(
            {recombination.Read().seq[i + recombination.DJInsertion().StartPosition() - 1],
             recombination.Read().seq[i + recombination.DJInsertion().StartPosition()]}));

    /** Total log probability */
    double result_probability = id_probs + del_probs + VD_non_genomic_ins_prob + DJ_non_genomic_ins_prob;
    return result_probability;
}

} // End namespace vdj_labeler
