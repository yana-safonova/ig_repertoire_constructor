#pragma once

#include "vj_alignment_info.hpp"
#include "vdj_hits.hpp"
#include "vdj_alignments/hits_calculator/d_hits_calculator/abstract_d_gene_hits_calculator.hpp"

namespace vdj_labeler {

class VDJHitsStorage {
private:
    std::vector <VDJHits> vdj_hits_;

public:
    VDJHitsStorage()                                 = default;
    VDJHitsStorage(const VDJHitsStorage&)            = default;
    VDJHitsStorage(VDJHitsStorage&&)                 = default;
    VDJHitsStorage& operator=(const VDJHitsStorage&) = default;
    VDJHitsStorage& operator=(VDJHitsStorage&&)      = default;

    explicit VDJHitsStorage(const vj_finder::VJAlignmentInfo &alignment_info) {
        for (const auto& alignment_record : alignment_info.AlignmentRecords()) {
            vdj_hits_.emplace_back(alignment_record);
        }
    }

    VDJHitsStorage(const vj_finder::VJAlignmentInfo &alignment_info,
                   AbstractDGeneHitsCalculator &d_gene_calculator)
    {
        for (const auto& alignment_record : alignment_info.AlignmentRecords()) {
            vdj_hits_.emplace_back(alignment_record, d_gene_calculator);
        }
    }

    void AddVDJHits(VDJHits vdj_hits) {
        vdj_hits_.emplace_back(std::move(vdj_hits));
    }

    size_t size() const { return vdj_hits_.size(); }

    typedef std::vector<VDJHits>::iterator vdj_hits_iterator;
    typedef std::vector<VDJHits>::const_iterator vdj_hits_const_iterator;

    vdj_hits_iterator       begin ()       { return vdj_hits_.begin (); }
    vdj_hits_const_iterator begin () const { return vdj_hits_.begin (); }
    vdj_hits_const_iterator cbegin() const { return vdj_hits_.cbegin(); }
    vdj_hits_iterator       end   ()       { return vdj_hits_.end   (); }
    vdj_hits_const_iterator end   () const { return vdj_hits_.end   (); }
    vdj_hits_const_iterator cend  () const { return vdj_hits_.cend  (); }

    // Users should not change the inside
    const VDJHits& operator[](const size_t &index) const {
        assert(index < size());
        return vdj_hits_[index];
    }
};

typedef std::shared_ptr <VDJHitsStorage> VDJHitsStoragePtr;

} // End namespace vdj_labeler
