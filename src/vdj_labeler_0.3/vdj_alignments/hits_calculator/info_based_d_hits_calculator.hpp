#pragma once

#include "vj_alignment_info.hpp"
#include "read_archive.hpp"
#include "germline_utils/germline_databases/immune_gene_database.hpp"
#include "../aligners/gene_segment_aligner.hpp"
#include "alignment_quality_checkers/alignment_quality_checker.hpp"
#include "../vdj_hits.hpp"

namespace vdj_labeler {

class InfoBasedDHitsCalculator {
    const core::ReadArchive &read_archive_;
    const vj_finder::VJAlignmentInfo &vj_alignment_info_;
    const germline_utils::ImmuneGeneDatabase &d_gene_database_;
    GeneSegmentAligner &d_gene_aligner_;
    vdj_labeler::AlignmentQualityChecker &quality_checker_;

    alignment_utils::AlignmentPositions ComputeDPositions(alignment_utils::ImmuneGeneAlignmentPositions v_positions,
                                                          alignment_utils::ImmuneGeneAlignmentPositions j_positions);

    alignment_utils::ImmuneGeneAlignmentPositions CreateDAlignmentPositions(
        alignment_utils::AlignmentPositions d_alignment_positions,
        germline_utils::ImmuneGenePtr gene_ptr,
        core::ReadPtr read_ptr);

    // possibly this methods should be replaced in a separate class
    bool DAlignmentPositionsAreGood(alignment_utils::AlignmentPositions d_alignment_positions);

public:
    InfoBasedDHitsCalculator(const core::ReadArchive &read_archive,
                             const vj_finder::VJAlignmentInfo &vj_alignment_info,
                             const germline_utils::ImmuneGeneDatabase &d_gene_database,
                             GeneSegmentAligner &d_gene_aligner,
                             vdj_labeler::AlignmentQualityChecker &quality_checker) :
        read_archive_(read_archive),
        vj_alignment_info_(vj_alignment_info),
        d_gene_database_(d_gene_database),
        d_gene_aligner_(d_gene_aligner),
        quality_checker_(quality_checker) { }

    // TODO Andrey: Implement, when I will see the usage.
    vdj_labeler::ImmuneGeneSegmentHitsPtr ComputeHits(core::ReadPtr read_ptr);
};

} // End namespace vdj_labeler