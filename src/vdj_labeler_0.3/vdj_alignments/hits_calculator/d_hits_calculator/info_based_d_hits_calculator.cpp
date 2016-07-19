#include "logger/logger.hpp"
#include "info_based_d_hits_calculator.hpp"

using namespace vdj_labeler;

alignment_utils::AlignmentPositions InfoBasedDHitsCalculator::ComputeDPositions(
    const ImmuneGeneSegmentHits &v_hits,
    const ImmuneGeneSegmentHits &j_hits) const
{
    size_t right_v_gene_pos = std::numeric_limits<size_t>::max();
    size_t left_g_gene_pos = 0;

    for (auto& v_hit : v_hits) {
        right_v_gene_pos = std::min(v_hit->EndQueryPosition(), right_v_gene_pos);
    }

    for (auto& j_hit : j_hits) {
        left_g_gene_pos = std::max(j_hit->StartQueryPosition(), left_g_gene_pos);
    }
    return alignment_utils::AlignmentPositions(
        std::make_pair<size_t, size_t>(right_v_gene_pos + 1, left_g_gene_pos - 1),
        std::make_pair<size_t, size_t>(0, size_t(-1)));
}

alignment_utils::ImmuneGeneAlignmentPositions InfoBasedDHitsCalculator::CreateDAlignmentPositions(
        alignment_utils::AlignmentPositions d_alignment_positions,
        const germline_utils::ImmuneGene& gene,
        const core::ReadPtr read_ptr) const
{
    d_alignment_positions.subject_pos.second = seqan::length(gene.seq()) - 1;
    assert(read_ptr != nullptr);
    return alignment_utils::ImmuneGeneAlignmentPositions(d_alignment_positions, gene, *read_ptr);
}

ImmuneGeneSegmentHitsPtr InfoBasedDHitsCalculator::ComputeDHits(const core::ReadPtr read_ptr,
                                                                const ImmuneGeneSegmentHits &v_hits,
                                                                const ImmuneGeneSegmentHits &j_hits) const
{
    auto d_positions = ComputeDPositions(v_hits, j_hits);

    ImmuneGeneSegmentHitsPtr d_hits_ptr(new ImmuneGeneSegmentHits(germline_utils::SegmentType::DiversitySegment,
                                                                  read_ptr));
    if (!d_alignment_position_checker_.DAlignmentPositionsAreGood(d_positions)) {
        TRACE("D positions are too short to generate alignment");
        TRACE(d_positions);
        // add single empty alignment and return hits storage
        seqan::Align<seqan::Dna5String> align;
        seqan::resize(seqan::rows(align), 2);
        d_hits_ptr->AddHit(alignment_utils::ImmuneGeneReadAlignmentPtr(new alignment_utils::ImmuneGeneReadAlignment(
            germline_utils::ImmuneGene(),
            *read_ptr,
            align)));
        return d_hits_ptr;
    }
    for (auto d_gene: d_gene_database_) {
        alignment_utils::ImmuneGeneAlignmentPositions d_alignment_pos = CreateDAlignmentPositions(
            d_positions,
            d_gene,
            read_ptr);
        auto d_alignment = d_gene_aligner_.ComputeAlignment(d_alignment_pos);
        if (quality_checker_.AlignmentIsGood(d_alignment)) {
            INFO(d_alignment->Alignment());
            d_hits_ptr->AddHit(d_alignment);
        }
    }
    return d_hits_ptr;
}
