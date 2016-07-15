#pragma once

#include "vj_alignment_info.hpp"
#include "read_archive.hpp"
#include "germline_utils/germline_databases/immune_gene_database.hpp"
#include "../aligners/gene_segment_aligner.hpp"
#include "alignment_quality_checkers/alignment_quality_checker.hpp"
#include "../vdj_hits.hpp"
#include "abstract_d_gene_hits_calculator.hpp"
#include "vdj_alignments/vdj_hits_storage.hpp"

namespace vdj_labeler {

class InfoBasedDHitsCalculator : AbstractDGeneHitsCalculator {
protected:
    VDJHitsStorage &vdj_hits_storage_;
    GeneSegmentAligner &d_gene_aligner_;

protected:
    alignment_utils::AlignmentPositions ComputeDPositions(size_t right_v, size_t left_j);

    alignment_utils::ImmuneGeneAlignmentPositions CreateDAlignmentPositions(
        alignment_utils::AlignmentPositions d_alignment_positions,
        germline_utils::ImmuneGenePtr gene_ptr,
        core::ReadPtr read_ptr);

    // possibly this methods should be replaced in a separate class
    bool DAlignmentPositionsAreGood(alignment_utils::AlignmentPositions d_alignment_positions);

public:
    InfoBasedDHitsCalculator(const core::ReadArchive &read_archive,
                             VDJHitsStorage &vdj_hits_storage,
                             const germline_utils::ImmuneGeneDatabase &d_gene_database,
                             GeneSegmentAligner &d_gene_aligner,
                             vdj_labeler::AlignmentQualityChecker &quality_checker) :
        AbstractDGeneHitsCalculator(read_archive, d_gene_database, quality_checker),
        vdj_hits_storage_(vdj_hits_storage),
        d_gene_aligner_(d_gene_aligner)
    { }

    // TODO Andrey: Implement, when I will see theHits(core::ReadPtr read_ptr);
    ImmuneGeneSegmentHitsPtr ComputeHits(core::ReadPtr read_ptr);
};

} // End namespace vdj_labeler