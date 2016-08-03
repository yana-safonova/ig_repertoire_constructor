#pragma once

#include "vdj_alignments/hits_calculator/alignment_quality_checkers/alignment_quality_checker.hpp"
#include "vj_alignment_structs.hpp"

namespace vdj_labeler {

// TODO Andrey: forward declaration here maybe is not ok, but I don't see any easy win here
class ImmuneGeneSegmentHits;
typedef std::shared_ptr<ImmuneGeneSegmentHits> ImmuneGeneSegmentHitsPtr;

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

    virtual ImmuneGeneSegmentHits ComputeDHits(const core::Read* read_ptr,
                                               const std::vector<vj_finder::VGeneHit> &v_hits,
                                               const std::vector<vj_finder::JGeneHit> &j_hits) const = 0;

    virtual ~AbstractDGeneHitsCalculator() { }
};

} // End namespace vdj_labeler