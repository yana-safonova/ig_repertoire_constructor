#include "compressed_cdr_set.hpp"

namespace cdr_labeler {
    void CompressedCDRSet::Initialize() {
        for(auto clone_it = clone_set_.cbegin(); clone_it != clone_set_.cend(); clone_it++) {
            if(clone_it->RegionIsEmpty(region_))
                continue;
            auto vj_hits = alignment_info_.GetVJHitsByRead(clone_it->Read());
            size_t index = compressed_cdrs_.size();
            CDRKey cdr_key(vj_hits.GetVHitByIndex(0).ImmuneGene().name(),
                           vj_hits.GetJHitByIndex(0).ImmuneGene().name(),
                           clone_it->GetRegionString(region_), index);
            if(compressed_cdrs_map_.find(cdr_key) == compressed_cdrs_map_.end()) {
                compressed_cdrs_.push_back(std::make_pair(cdr_key, 1));
                compressed_cdrs_map_[cdr_key] =  index;
            }
            else {
                index = compressed_cdrs_map_.at(cdr_key);
                compressed_cdrs_[index].second++;
            }
        }
        size_t max_abundance = 0;
        for(auto it = cbegin(); it != cend(); it++)
            max_abundance = std::max<size_t>(max_abundance, it->second);
        INFO(compressed_cdrs_.size() << " unique " << region_ << " sequences were created, max abundance: " <<
                     max_abundance);
    }
}