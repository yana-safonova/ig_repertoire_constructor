#pragma once

#include "vj_alignment_info.hpp"
#include "vdj_hits.hpp"

namespace vdj_labeler {

class VDJHitsStorage {
    std::vector <VDJHitsPtr> vdj_hits_;

public:

    VDJHitsStorage(const vj_finder::VJAlignmentInfo &alignment_info)
    {
        for (auto& alignment_record : alignment_info.AlignmentRecords()) {
            vdj_hits_.emplace_back(std::make_shared<VDJHits>(VDJHits(alignment_record)));
        }
    }

    void AddVDJHits(const VDJHitsPtr &vdj_hits_ptr) {
        vdj_hits_.emplace_back(vdj_hits_ptr);
    }

    size_t size() const { return vdj_hits_.size(); }

    typedef std::vector<VDJHitsPtr>::const_iterator vdj_hits_citerator;

    vdj_hits_citerator cbegin() const { return vdj_hits_.cbegin(); }

    vdj_hits_citerator cend() const { return vdj_hits_.cend(); }

    VDJHitsPtr operator[](const size_t &index) const {
        assert(index < size());
        return vdj_hits_[index];
    }
};

typedef std::shared_ptr <VDJHitsStorage> VDJHitsStoragePtr;

} // End namespace vdj_labeler
