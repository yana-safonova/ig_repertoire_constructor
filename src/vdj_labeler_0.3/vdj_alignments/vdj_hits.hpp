#pragma once

#include <memory>

#include "vj_alignment_structs.hpp"
#include "alignment_utils/alignment_positions.hpp"
#include "alignment_utils/pairwise_alignment.hpp"
#include "germline_utils/germline_gene_type.hpp"
#include "read_archive.hpp"
#include "vdj_alignments/hits_calculator/d_hits_calculator/abstract_d_gene_hits_calculator.hpp"

namespace vdj_labeler {

class ImmuneGeneSegmentHits {
protected:
    germline_utils::SegmentType segment_type_;
    const core::Read* read_ptr_;
    std::vector<alignment_utils::ImmuneGeneReadAlignment> hits_;

public:
    ImmuneGeneSegmentHits(const germline_utils::SegmentType &gene_type = germline_utils::SegmentType(),
                          const core::Read* read_ptr = nullptr);
    ImmuneGeneSegmentHits(const germline_utils::SegmentType &gene_type,
                          const core::Read &read);
    ImmuneGeneSegmentHits(const germline_utils::SegmentType &gene_type,
                          const core::Read* read_ptr,
                          const std::vector<vj_finder::ImmuneGeneHit>&);

    ImmuneGeneSegmentHits(const ImmuneGeneSegmentHits&)            = default;
    ImmuneGeneSegmentHits(ImmuneGeneSegmentHits&&)                 = default;
    ImmuneGeneSegmentHits& operator=(const ImmuneGeneSegmentHits&) = default;
    ImmuneGeneSegmentHits& operator=(ImmuneGeneSegmentHits&&)      = default;

    void AddHit(alignment_utils::ImmuneGeneReadAlignment hit);

    size_t size() const;

    typedef std::vector<alignment_utils::ImmuneGeneReadAlignment>::iterator hits_iterator;
    typedef std::vector<alignment_utils::ImmuneGeneReadAlignment>::const_iterator hits_citerator;

    hits_iterator  begin ()       { return hits_.begin (); }
    hits_citerator begin () const { return hits_.begin (); }
    hits_citerator cbegin() const { return hits_.cbegin(); }
    hits_iterator  end   ()       { return hits_.end   (); }
    hits_citerator end   () const { return hits_.end   (); }
    hits_citerator cend  () const { return hits_.cend  (); }

    const alignment_utils::ImmuneGeneReadAlignment& operator[](const size_t &index) const;

    const germline_utils::SegmentType& GeneType() const;

    const core::Read* ReadPtr() const { return read_ptr_; }
};

typedef std::shared_ptr<ImmuneGeneSegmentHits> ImmuneGeneSegmentHitsPtr;

//------------------------------------------------------------

class VGeneHits : public ImmuneGeneSegmentHits {
public:
    VGeneHits(const germline_utils::SegmentType &gene_type = germline_utils::SegmentType(),
          const core::Read* read_ptr = nullptr) :
        ImmuneGeneSegmentHits(gene_type, read_ptr)
    { }

    VGeneHits(const germline_utils::SegmentType &gene_type,
          const core::Read &read) :
        ImmuneGeneSegmentHits(gene_type, read)
    { }

    VGeneHits(const germline_utils::SegmentType &gene_type,
          const core::Read* read_ptr,
          const std::vector<vj_finder::VGeneHit>& hits);

    VGeneHits(const VGeneHits&)            = default;
    VGeneHits(VGeneHits&&)                 = default;
    VGeneHits& operator=(const VGeneHits&) = default;
    VGeneHits& operator=(VGeneHits&&)      = default;
};

//------------------------------------------------------------

class JGeneHits : public ImmuneGeneSegmentHits {
public:
    JGeneHits(const germline_utils::SegmentType &gene_type = germline_utils::SegmentType(),
          const core::Read* read_ptr = nullptr) :
        ImmuneGeneSegmentHits(gene_type, read_ptr)
    { }

    JGeneHits(const germline_utils::SegmentType &gene_type,
          const core::Read &read) :
        ImmuneGeneSegmentHits(gene_type, read)
    { }

    JGeneHits(const germline_utils::SegmentType &gene_type,
          const core::Read* read_ptr,
          const std::vector<vj_finder::JGeneHit>& hits);

    JGeneHits(const JGeneHits&)            = default;
    JGeneHits(JGeneHits&&)                 = default;
    JGeneHits& operator=(const JGeneHits&) = default;
    JGeneHits& operator=(JGeneHits&&)      = default;
};

//------------------------------------------------------------

typedef ImmuneGeneSegmentHits DGeneHits;

//------------------------------------------------------------

class VDJHits {
private:
    const core::Read* read_ptr_;
    VGeneHits v_hits_;
    DGeneHits d_hits_;
    JGeneHits j_hits_;

public:
    VDJHits(const core::Read* read_ptr = nullptr);

    VDJHits(const core::Read* read_ptr,
            const std::vector<vj_finder::VGeneHit> &v_hits,
            const std::vector<vj_finder::JGeneHit> &j_hits);

    VDJHits(const core::Read* read_ptr,
            const std::vector<vj_finder::VGeneHit> &v_hits,
            const std::vector<vj_finder::JGeneHit> &j_hits,
            AbstractDGeneHitsCalculator &d_gene_calculator);

    VDJHits(const vj_finder::VJHits &vj_hits);

    VDJHits(const vj_finder::VJHits &vj_hits, AbstractDGeneHitsCalculator &d_gene_calculator);

    VDJHits(const VDJHits&)            = default;
    VDJHits(VDJHits&&)                 = default;
    VDJHits& operator=(const VDJHits&) = default;
    VDJHits& operator=(VDJHits&&)      = default;

    void AddIgGeneAlignment(const germline_utils::SegmentType &gene_type,
                            alignment_utils::ImmuneGeneReadAlignment alignment);

    void AddIgGeneAlignment(alignment_utils::ImmuneGeneReadAlignment alignment);

    size_t VHitsNumber() const { return v_hits_.size(); }
    size_t DHitsNumber() const { return d_hits_.size(); }
    size_t JHitsNumber() const { return j_hits_.size(); }

    const VGeneHits &VHits() const { return v_hits_; }
    const DGeneHits &DHits() const { return d_hits_; }
    const JGeneHits &JHits() const { return j_hits_; }

    void SetVHits(const VGeneHits &v_hits) { v_hits_ = v_hits; }
    void SetDHits(const DGeneHits &d_hits) { d_hits_ = d_hits; }
    void SetJHits(const JGeneHits &j_hits) { j_hits_ = j_hits; }

    const core::Read* ReadPtr() const { return read_ptr_; }
    const core::Read& Read() const {
        assert(read_ptr_ != nullptr);
        return *read_ptr_;
    }

    const alignment_utils::ImmuneGeneReadAlignment& GetAlignmentByIndex(const germline_utils::SegmentType &gene_type,
                                                                        const size_t &index);
};

typedef std::shared_ptr<VDJHits> VDJHitsPtr;

} // End namespace vdj_labeler
