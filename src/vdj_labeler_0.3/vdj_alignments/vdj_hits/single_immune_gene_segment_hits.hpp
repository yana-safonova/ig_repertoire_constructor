//
// Created by Andrew Bzikadze on 8/5/16.
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

class SingleImmuneGeneSegmentHits {
private:
    germline_utils::SegmentType segment_type_;
    const core::Read* read_ptr_;
    std::vector<alignment_utils::ImmuneGeneReadAlignment> hits_;

public:
    SingleImmuneGeneSegmentHits(const germline_utils::SegmentType &gene_type = germline_utils::SegmentType(),
                                const core::Read* read_ptr = nullptr) :
        segment_type_(gene_type),
        read_ptr_(read_ptr)
    { }

    SingleImmuneGeneSegmentHits(const germline_utils::SegmentType &gene_type,
                                const core::Read &read) :
        SingleImmuneGeneSegmentHits(gene_type, &read)
    { }

    // vj_finder_ImmuneGeneHit could be vj::finder {ImmuneGeneHit, VGeneHit, JGeneHit}.
    template<class vj_finder_ImmuneGeneHit>
    SingleImmuneGeneSegmentHits(const germline_utils::SegmentType &segment_type,
                                const core::Read* read_ptr,
                                const std::vector<vj_finder_ImmuneGeneHit>& hits) :
        SingleImmuneGeneSegmentHits(segment_type, read_ptr)
    {
        vj_finder::ImmuneGeneAlignmentConverter converter;
        assert(read_ptr != nullptr);
        for (const auto& hit : hits) {
            hits_.emplace_back(converter.ConvertToAlignment(hit.ImmuneGene(), hit.Read(), hit.BlockAlignment()));
        }
    }

    SingleImmuneGeneSegmentHits(const SingleImmuneGeneSegmentHits&)            = default;
    SingleImmuneGeneSegmentHits(SingleImmuneGeneSegmentHits&&)                 = default;
    SingleImmuneGeneSegmentHits& operator=(const SingleImmuneGeneSegmentHits&) = default;
    SingleImmuneGeneSegmentHits& operator=(SingleImmuneGeneSegmentHits&&)      = default;

    void AddHit(alignment_utils::ImmuneGeneReadAlignment hit);

    size_t size() const { return hits_.size(); }

    typedef std::vector<alignment_utils::ImmuneGeneReadAlignment>::iterator hits_iterator;
    typedef std::vector<alignment_utils::ImmuneGeneReadAlignment>::const_iterator hits_citerator;

    hits_iterator  begin ()       { return hits_.begin (); }
    hits_citerator begin () const { return hits_.begin (); }
    hits_citerator cbegin() const { return hits_.cbegin(); }
    hits_iterator  end   ()       { return hits_.end   (); }
    hits_citerator end   () const { return hits_.end   (); }
    hits_citerator cend  () const { return hits_.cend  (); }

    const alignment_utils::ImmuneGeneReadAlignment& operator[](const size_t &index) const;

    const germline_utils::SegmentType& GeneType() const { return segment_type_; }

    const core::Read& Read() const {
        assert(read_ptr_ != nullptr);
        return *read_ptr_;
    }

    const core::Read* ReadPtr() const { return read_ptr_; }
};

} // End namespace vdj_labeler