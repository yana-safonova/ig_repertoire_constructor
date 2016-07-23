#include "logger/logger.hpp"
#include "info_based_d_hits_calculator.hpp"

using namespace vdj_labeler;

alignment_utils::AlignmentPositions InfoBasedDHitsCalculator::ComputeDPositions(
    const ImmuneGeneSegmentHits &v_hits,
    const ImmuneGeneSegmentHits &j_hits) const
{
    size_t right_v_gene_pos = std::numeric_limits<size_t>::max();
    size_t left_g_gene_pos = 0;

    for (const auto& v_hit : v_hits) {
        right_v_gene_pos = std::min(v_hit.EndQueryPosition(), right_v_gene_pos);
    }

    for (const auto& j_hit : j_hits) {
        left_g_gene_pos = std::max(j_hit.StartQueryPosition(), left_g_gene_pos);
    }
    return alignment_utils::AlignmentPositions(
        std::make_pair<size_t, size_t>(right_v_gene_pos + 1, left_g_gene_pos - 1),
        std::make_pair<size_t, size_t>(0, size_t(-1)));
}

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
                                                             const ImmuneGeneSegmentHits &v_hits,
                                                             const ImmuneGeneSegmentHits &j_hits) const
{
    auto d_positions = ComputeDPositions(v_hits, j_hits);

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
    // INFO("\n ind for");
    // for (size_t i = 0; i < d_hits.size(); ++i) {
    //     INFO((d_hits[i]->subject().Segment() == germline_utils::SegmentType::DiversitySegment));
    // }
    // for (size_t i = 0; i < d_hits.size(); ++i) {
    //     INFO(d_hits[i]->subject().seq());
    // }
    // INFO("\n ord for");
    // for (auto it = d_hits.begin(); it != d_hits.end(); ++it) {
    //     INFO(((*it)->subject().Segment() == germline_utils::SegmentType::DiversitySegment));
    // }
    // INFO("\n range-based for");
    // for (const auto& d_hit : d_hits) {
    //     INFO((d_hit->subject().Segment() == germline_utils::SegmentType::DiversitySegment));
    // }
    // INFO("\n");
    return d_hits;
}
