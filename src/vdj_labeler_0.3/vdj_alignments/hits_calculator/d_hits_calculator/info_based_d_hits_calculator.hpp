#pragma once

#include "vj_alignment_info.hpp"
#include "read_archive.hpp"
#include "germline_utils/germline_databases/immune_gene_database.hpp"
#include "vdj_alignments/aligners/gene_segment_aligner.hpp"
#include "vdj_alignments/hits_calculator/alignment_quality_checkers/alignment_quality_checker.hpp"
#include "vdj_alignments/vdj_hits.hpp"
#include "abstract_d_gene_hits_calculator.hpp"
#include "vdj_alignments/vdj_hits_storage.hpp"
#include "vdj_alignments/hits_calculator/d_alignment_positions_checkers/abstract_d_alignment_position_checker.hpp"
#include "vdj_alignments/hits_calculator/d_alignment_positions_calculator/abstract_d_alignment_positions_calculator.hpp"

namespace vdj_labeler {

class InfoBasedDHitsCalculator : public AbstractDGeneHitsCalculator {
protected:
    GeneSegmentAligner &d_gene_aligner_;
    AbstractDAlignmentPositionChecker &d_alignment_position_checker_;
    AbstractDAlignmentPositionsCalculator &d_alignment_positions_calculator_;

protected:
    alignment_utils::ImmuneGeneAlignmentPositions CreateDAlignmentPositions(
        alignment_utils::AlignmentPositions d_alignment_positions,
        const germline_utils::ImmuneGene* gene_ptr,
        const core::Read* read_ptr) const;

public:
    InfoBasedDHitsCalculator(const germline_utils::ImmuneGeneDatabase &d_gene_database,
                             GeneSegmentAligner &d_gene_aligner,
                             AlignmentQualityChecker &quality_checker,
                             AbstractDAlignmentPositionChecker &d_alignment_position_checker,
                             AbstractDAlignmentPositionsCalculator &d_alignment_position_calculator) :
        AbstractDGeneHitsCalculator(d_gene_database, quality_checker),
        d_gene_aligner_(d_gene_aligner),
        d_alignment_position_checker_(d_alignment_position_checker),
        d_alignment_positions_calculator_(d_alignment_position_calculator)
    { }

    virtual ImmuneGeneSegmentHits ComputeDHits(const core::Read* read_ptr,
                                               const std::vector<vj_finder::VGeneHit> &v_hits,
                                               const std::vector<vj_finder::JGeneHit> &j_hits) const override;
};

} // End namespace vdj_labeler