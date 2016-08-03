//
// Created by Andrew Bzikadze on 7/26/16.
//

#include "custom_d_alignment_positions_calculator.hpp"

using namespace vdj_labeler;

alignment_utils::AlignmentPositions CustomDAlignmentPositionsCalculator::ComputeDAlignmentPositions(
        const std::vector<vj_finder::VGeneHit> &v_hits,
        const std::vector<vj_finder::JGeneHit> &j_hits) const
{
    size_t left_j_gene_pos = std::numeric_limits<size_t>::max();
    size_t right_v_gene_pos = 0;

    for (const auto& v_hit : v_hits) {
        right_v_gene_pos = std::max(v_hit.BlockAlignment().last_match_read_pos(), right_v_gene_pos);
    }

    for (const auto& j_hit : j_hits) {
        left_j_gene_pos = std::min(j_hit.BlockAlignment().first_match_read_pos(), left_j_gene_pos);
    }
    return alignment_utils::AlignmentPositions(
        std::make_pair<size_t, size_t>(right_v_gene_pos + 1, left_j_gene_pos - 1),
        std::make_pair<size_t, size_t>(0, size_t(-1)));
}
