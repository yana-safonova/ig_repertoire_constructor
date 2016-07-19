#pragma once

#include <memory>

#include "vj_alignment_structs.hpp"
#include "alignment_utils/alignment_positions.hpp"
#include "alignment_utils/pairwise_alignment.hpp"
#include "germline_utils/germline_gene_type.hpp"
#include "read_archive.hpp"
#include "vdj_alignments/hits_calculator/d_hits_calculator/abstract_d_gene_hits_calculator.hpp"

namespace vdj_labeler {

typedef std::vector<alignment_utils::ImmuneGeneReadAlignmentPtr>::iterator hits_iterator;
typedef std::vector<alignment_utils::ImmuneGeneReadAlignmentPtr>::const_iterator hits_citerator;

class ImmuneGeneSegmentHits {
    germline_utils::SegmentType segment_type_;
    core::ReadPtr read_ptr_;
    std::vector<alignment_utils::ImmuneGeneReadAlignmentPtr> hits_;

public:
    ImmuneGeneSegmentHits(const germline_utils::SegmentType &gene_type, const core::ReadPtr &read_ptr);

    ImmuneGeneSegmentHits(const germline_utils::SegmentType &gene_type, const core::ReadPtr &read_ptr,
                          const std::vector<vj_finder::ImmuneGeneHitPtr>&);

    void AddHit(const alignment_utils::ImmuneGeneReadAlignmentPtr &hit);

    size_t size() const;

    hits_iterator  begin ();
    hits_citerator begin () const;
    hits_citerator cbegin() const;
    hits_iterator  end   ();
    hits_citerator end   () const;
    hits_citerator cend  () const;

    alignment_utils::ImmuneGeneReadAlignmentPtr operator[](const size_t &index);

    germline_utils::SegmentType GeneType() const;
};

typedef std::shared_ptr<ImmuneGeneSegmentHits> ImmuneGeneSegmentHitsPtr;

//------------------------------------------------------------

class VDJHits {
private:
    core::ReadPtr read_ptr_;
    ImmuneGeneSegmentHits v_hits_;
    ImmuneGeneSegmentHits d_hits_;
    ImmuneGeneSegmentHits j_hits_;

public:
    VDJHits(const core::ReadPtr &read_ptr);

    VDJHits(const core::ReadPtr &read_ptr,
            const std::vector<vj_finder::ImmuneGeneHitPtr>& v_hits,
            const std::vector<vj_finder::ImmuneGeneHitPtr>& j_hits);

    VDJHits(const core::ReadPtr &read_ptr,
            const std::vector<vj_finder::ImmuneGeneHitPtr>& v_hits,
            const std::vector<vj_finder::ImmuneGeneHitPtr>& j_hits,
            const AbstractDGeneHitsCalculator &d_gene_calculator);

    VDJHits(const vj_finder::VJHits &vj_hits);

    VDJHits(const vj_finder::VJHits &vj_hits, const AbstractDGeneHitsCalculator &d_gene_calculator);

    void AddIgGeneAlignment(const germline_utils::SegmentType &gene_type,
                            const alignment_utils::ImmuneGeneReadAlignmentPtr &alignment_ptr);

    void AddIgGeneAlignment(const alignment_utils::ImmuneGeneReadAlignmentPtr &alignment_ptr);

    size_t VHitsNumber() const { return v_hits_.size(); }

    size_t DHitsNumber() const { return d_hits_.size(); }

    size_t JHitsNumber() const { return j_hits_.size(); }

    const ImmuneGeneSegmentHits &VHits() const { return v_hits_; }
    const ImmuneGeneSegmentHits &DHits() const { return d_hits_; }
    const ImmuneGeneSegmentHits &JHits() const { return j_hits_; }

    void SetVHits(const ImmuneGeneSegmentHits &v_hits_) { VDJHits::v_hits_ = v_hits_; }
    void SetDHits(const ImmuneGeneSegmentHits &d_hits_) { VDJHits::d_hits_ = d_hits_; }
    void SetJHits(const ImmuneGeneSegmentHits &j_hits_) { VDJHits::j_hits_ = j_hits_; }

    core::ReadPtr Read() const { return read_ptr_; }

    alignment_utils::ImmuneGeneReadAlignmentPtr GetAlignmentByIndex(const germline_utils::SegmentType &gene_type,
                                                                    const size_t &index);
};

typedef std::shared_ptr<VDJHits> VDJHitsPtr;

} // End namespace vdj_labeler
