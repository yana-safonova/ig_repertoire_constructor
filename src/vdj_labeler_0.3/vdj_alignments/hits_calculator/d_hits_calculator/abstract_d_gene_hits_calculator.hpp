#pragma once

#include "vdj_alignments/hits_calculator/alignment_quality_checkers/alignment_quality_checker.hpp"
#include "vj_alignment_structs.hpp"
#include "vdj_alignments/vdj_hits/d_gene_segment_hit.hpp"

namespace vdj_labeler {

class AbstractDGeneHitsCalculator {
protected:
    const germline_utils::ImmuneGeneDatabase &d_gene_database_;
    AlignmentQualityChecker &quality_checker_;

public:
    AbstractDGeneHitsCalculator() = delete;
    AbstractDGeneHitsCalculator(const AbstractDGeneHitsCalculator &) = delete;
    AbstractDGeneHitsCalculator& operator=(const AbstractDGeneHitsCalculator &) = delete;
    AbstractDGeneHitsCalculator(AbstractDGeneHitsCalculator &&) = delete;
    AbstractDGeneHitsCalculator& operator=(AbstractDGeneHitsCalculator&&) = delete;

    AbstractDGeneHitsCalculator(const germline_utils::ImmuneGeneDatabase &d_gene_database,
                                AlignmentQualityChecker &quality_checker) :
        d_gene_database_(d_gene_database),
        quality_checker_(quality_checker)
    { }

    virtual DGeneHits ComputeDHits(const core::Read* read_ptr,
                                   const std::vector<vj_finder::VGeneHit> &v_hits,
                                   const std::vector<vj_finder::JGeneHit> &j_hits) const = 0;

    virtual ~AbstractDGeneHitsCalculator() { }
};

} // End namespace vdj_labeler