#include "logger/logger.hpp"
#include "info_based_d_hits_calculator.hpp"

using namespace vdj_labeler;

alignment_utils::ImmuneGeneAlignmentPositions InfoBasedDHitsCalculator::CreateDAlignmentPositions(
        alignment_utils::AlignmentPositions d_alignment_positions,
        const germline_utils::ImmuneGene* gene_ptr,
        const core::Read* read_ptr) const
{
    d_alignment_positions.subject_pos.second = seqan::length(gene_ptr->seq()) - 1;
    assert(read_ptr != nullptr);
    return alignment_utils::ImmuneGeneAlignmentPositions(std::move(d_alignment_positions), gene_ptr, read_ptr);
}

ImmuneGeneSegmentHits InfoBasedDHitsCalculator::ComputeDHits(const core::Read* read_ptr,
                                           const std::vector<vj_finder::VGeneHit> &v_hits,
                                           const std::vector<vj_finder::JGeneHit> &j_hits) const
{
    auto d_positions = d_alignment_positions_calculator_.ComputeDAlignmentPositions(v_hits, j_hits);

    ImmuneGeneSegmentHits d_hits(germline_utils::SegmentType::DiversitySegment, read_ptr);
    if (!d_alignment_position_checker_.DAlignmentPositionsAreGood(d_positions)) {
        TRACE("D positions are too short to generate alignment");
        TRACE(d_positions);
        // add single empty alignment and return hits storage
        seqan::Align<seqan::Dna5String> align;
        seqan::resize(seqan::rows(align), 2);
        d_hits.AddHit(alignment_utils::ImmuneGeneReadAlignment(nullptr, read_ptr, align, 0));
        return d_hits;
    }
    for (const auto& d_gene: d_gene_database_) {
        alignment_utils::ImmuneGeneAlignmentPositions d_alignment_pos =
            CreateDAlignmentPositions(d_positions, &d_gene, read_ptr);
        auto d_alignment = d_gene_aligner_.ComputeAlignment(d_alignment_pos);
        if (quality_checker_.AlignmentIsGood(d_alignment)) {
            // INFO(d_alignment->Alignment());
            d_hits.AddHit(std::move(d_alignment));
        }
    }
    return d_hits;
}
