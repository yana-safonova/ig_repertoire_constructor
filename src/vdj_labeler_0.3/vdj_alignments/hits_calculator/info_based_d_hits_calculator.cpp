#include "logger/logger.hpp"
#include "info_based_d_hits_calculator.hpp"
#include "vdj_alignments/vdj_hits.hpp"

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>

namespace vdj_labeler {

alignment_utils::AlignmentPositions InfoBasedDHitsCalculator::ComputeDPositions(size_t right_v, size_t left_j) {
    return alignment_utils::AlignmentPositions(std::make_pair<size_t, size_t>(right_v + 1, left_j - 1),
                                               std::make_pair<size_t, size_t>(0, size_t(-1)));
}

alignment_utils::ImmuneGeneAlignmentPositions InfoBasedDHitsCalculator::CreateDAlignmentPositions(
        alignment_utils::AlignmentPositions d_alignment_positions,
        germline_utils::ImmuneGenePtr gene_ptr,
        core::ReadPtr read_ptr)
{
    d_alignment_positions.subject_pos.second = seqan::length(gene_ptr->seq()) - 1;
    assert(gene_ptr != nullptr);
    assert(read_ptr != nullptr);
    return alignment_utils::ImmuneGeneAlignmentPositions(d_alignment_positions, *gene_ptr, *read_ptr);
}

// todo: refactor magic number!
bool InfoBasedDHitsCalculator::DAlignmentPositionsAreGood(alignment_utils::AlignmentPositions d_alignment_positions) {
    return d_alignment_positions.QueryAlignmentLength() >= 5;
}

vdj_labeler::ImmuneGeneSegmentHitsPtr InfoBasedDHitsCalculator::ComputeHits(core::ReadPtr read_ptr) {
    size_t read_index = read_archive_.GetIndexByReadName(read_ptr->name);
    assert(read_index == read_ptr->id);

    const VDJHitsPtr& vdj_hits_ptr = vdj_hits_storage_[read_index];

    size_t right_v_gene_pos = std::numeric_limits<size_t>::max();
    size_t left_g_gene_pos = 0;

    for (auto& v_hit : vdj_hits_ptr->VHits()) {
        right_v_gene_pos = std::min(v_hit->EndQueryPosition(), right_v_gene_pos);
    }

    for (auto& j_hit : vdj_hits_ptr->JHits()) {
        left_g_gene_pos = std::max(j_hit->StartQueryPosition(), left_g_gene_pos);
    }
    auto d_positions = ComputeDPositions(right_v_gene_pos, left_g_gene_pos);

    ImmuneGeneSegmentHitsPtr d_hits_ptr(new vdj_labeler::ImmuneGeneSegmentHits(
                                                                        germline_utils::SegmentType::DiversitySegment,
                                                                        read_ptr));
    if (!DAlignmentPositionsAreGood(d_positions)) {
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
            std::make_shared<germline_utils::ImmuneGene>(d_gene),
            read_ptr);
        auto d_alignment = d_gene_aligner_.ComputeAlignment(d_alignment_pos);
        if (quality_checker_.AlignmentIsGood(d_alignment))
            d_hits_ptr->AddHit(d_alignment);
    }
    return d_hits_ptr;
}

} // End namespace vdj_labeler