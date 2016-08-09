//
// Created by Andrew Bzikadze on 8/6/16.
//

#pragma once

#include "vj_alignment_structs.hpp"
#include "alignment_utils/alignment_positions.hpp"
#include "alignment_utils/pairwise_alignment.hpp"
#include "germline_utils/germline_gene_type.hpp"
#include "read_archive.hpp"
#include "immune_gene_alignment_converter.hpp"

namespace vdj_labeler {

// DGeneHit is a vector of alignments of D genes. This class covers tandem D genes.
class DGeneHit {
private:
    const core::Read* read_ptr_;
    std::vector<alignment_utils::ImmuneGeneReadAlignment> d_genes_hits_;

public:
    DGeneHit(const core::Read* read_ptr = nullptr) :
        read_ptr_(read_ptr)
    { }

    DGeneHit(const core::Read &read) :
        DGeneHit(&read)
    { }

    DGeneHit(const core::Read* read_ptr,
             std::vector<alignment_utils::ImmuneGeneReadAlignment> d_genes_hits) :
        read_ptr_(read_ptr),
        d_genes_hits_(std::move(d_genes_hits))
    { }

    DGeneHit(const DGeneHit&)            = default;
    DGeneHit(DGeneHit&&)                 = default;
    DGeneHit& operator=(const DGeneHit&) = default;
    DGeneHit& operator=(DGeneHit&&)      = default;

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
        VERIFY(index < size());
        return d_genes_hits_[index];
    }

    const core::Read& Read() const {
        VERIFY(read_ptr_ != nullptr);
        return *read_ptr_;
    }

    const core::Read* ReadPtr() const { return read_ptr_; }
};

//--------------------------------------------------------------------------------/

// DGeneHits is a vector of DGeneHit-s. This class is a D gene -analog of SingleImmuneGeneSegment class.
class DGeneHits {
private:
    const core::Read* read_ptr_;
    std::vector<DGeneHit> hits_;

public:
    DGeneHits(const core::Read* read_ptr = nullptr) :
        read_ptr_(read_ptr)
    { }

    DGeneHits(const core::Read &read) :
        DGeneHits(&read)
    { }

    DGeneHits(const DGeneHits&)            = default;
    DGeneHits(DGeneHits&&)                 = default;
    DGeneHits& operator=(const DGeneHits&) = default;
    DGeneHits& operator=(DGeneHits&&)      = default;

    size_t size() const { return hits_.size(); }

    typedef std::vector<DGeneHit>::iterator hits_iterator;
    typedef std::vector<DGeneHit>::const_iterator hits_citerator;

    hits_iterator  begin ()       { return hits_.begin (); }
    hits_citerator begin () const { return hits_.begin (); }
    hits_citerator cbegin() const { return hits_.cbegin(); }
    hits_iterator  end   ()       { return hits_.end   (); }
    hits_citerator end   () const { return hits_.end   (); }
    hits_citerator cend  () const { return hits_.cend  (); }

    const DGeneHit& operator[](const size_t &index) const {
        VERIFY(index < size());
        return hits_[index];
    }

    const core::Read& Read() const {
        VERIFY(read_ptr_ != nullptr);
        return *read_ptr_;
    }

    const core::Read* ReadPtr() const { return read_ptr_; }

    void AddHit(DGeneHit hit) {
        hits_.emplace_back(std::move(hit));
    }
};

} // End vdj_labeler
