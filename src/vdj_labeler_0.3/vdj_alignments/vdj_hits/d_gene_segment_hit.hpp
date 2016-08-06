//
// Created by Andrew Bzikadze on 8/6/16.
//

#pragma once

#include "vj_alignment_structs.hpp"
#include "alignment_utils/alignment_positions.hpp"
#include "alignment_utils/pairwise_alignment.hpp"
#include "germline_utils/germline_gene_type.hpp"
#include "read_archive.hpp"
#include "vdj_alignments/hits_calculator/d_hits_calculator/abstract_d_gene_hits_calculator.hpp"
#include "immune_gene_alignment_converter.hpp"

namespace vdj_labeler {

class DGeneSegmentHit {
private:
    const core::Read* read_ptr_;
    std::vector<alignment_utils::ImmuneGeneReadAlignment> d_genes_hits_;

public:
    DGeneSegmentHit(const core::Read* read_ptr = nullptr) :
        read_ptr_(read_ptr)
    { }

    DGeneSegmentHit(const core::Read &read) :
        DGeneSegmentHit(&read)
    { }

    DGeneSegmentHit(const DGeneSegmentHit&)            = default;
    DGeneSegmentHit(DGeneSegmentHit&&)                 = default;
    DGeneSegmentHit& operator=(const DGeneSegmentHit&) = default;
    DGeneSegmentHit& operator=(DGeneSegmentHit&&)      = default;

    size_t size() const { return d_genes_hits_.size(); }

    typedef std::vector<alignment_utils::ImmuneGeneReadAlignment>::iterator hits_iterator;
    typedef std::vector<alignment_utils::ImmuneGeneReadAlignment>::const_iterator hits_citerator;

    hits_iterator  begin ()       { return d_genes_hits_.begin (); }
    hits_citerator begin () const { return d_genes_hits_.begin (); }
    hits_citerator cbegin() const { return d_genes_hits_.cbegin(); }
    hits_iterator  end   ()       { return d_genes_hits_.end   (); }
    hits_citerator end   () const { return d_genes_hits_.end   (); }
    hits_citerator cend  () const { return d_genes_hits_.cend  (); }

    const alignment_utils::ImmuneGeneReadAlignment& operator[](const size_t &index) const {
        assert(index < size());
        return d_genes_hits_[index];
    }

    const core::Read& Read() const {
        assert(read_ptr_ != nullptr);
        return *read_ptr_;
    }

    const core::Read* ReadPtr() const { return read_ptr_; }
};

//--------------------------------------------------------------------------------/

class DGeneSegmentHits {
private:
    const core::Read* read_ptr_;
    std::vector<DGeneSegmentHit> hits_;

public:
    DGeneSegmentHits(const core::Read* read_ptr = nullptr) :
        read_ptr_(read_ptr)
    { }

    DGeneSegmentHits(const core::Read &read) :
        DGeneSegmentHits(&read)
    { }

    DGeneSegmentHits(const DGeneSegmentHits&)            = default;
    DGeneSegmentHits(DGeneSegmentHits&&)                 = default;
    DGeneSegmentHits& operator=(const DGeneSegmentHits&) = default;
    DGeneSegmentHits& operator=(DGeneSegmentHits&&)      = default;

    size_t size() const { return hits_.size(); }

    typedef std::vector<DGeneSegmentHit>::iterator hits_iterator;
    typedef std::vector<DGeneSegmentHit>::const_iterator hits_citerator;

    hits_iterator  begin ()       { return hits_.begin (); }
    hits_citerator begin () const { return hits_.begin (); }
    hits_citerator cbegin() const { return hits_.cbegin(); }
    hits_iterator  end   ()       { return hits_.end   (); }
    hits_citerator end   () const { return hits_.end   (); }
    hits_citerator cend  () const { return hits_.cend  (); }

    const DGeneSegmentHit& operator[](const size_t &index) const {
        assert(index < size());
        return hits_[index];
    }

    const core::Read& Read() const {
        assert(read_ptr_ != nullptr);
        return *read_ptr_;
    }

    const core::Read* ReadPtr() const { return read_ptr_; }
};

} // End vdj_labeler
